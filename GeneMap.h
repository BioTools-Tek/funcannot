#ifndef GENEMAP_H
#define GENEMAP_H

#include "ProteinHandler.h"

typedef uint t_exon;
//typedef QList<t_exon> exondata;

//FOR ALL
static uint countlines(QString file){
    uint lines=0;
    //For lines in infile
    QFile inputFile(file);
    if (inputFile.open(QIODevice::ReadOnly))
    {
        while (!inputFile.atEnd()){
            inputFile.readLine();
            lines ++;
        }
    }
    return lines;
    inputFile.close();
}

static void progress(const char * filename, int tab, int progress){
    cerr << '\r' << flush;
    for (int tn=0; tn < (2*tab); tn++){
        cerr << '\t' << flush;
    }
    cerr << '#' << filename << ':' << progress << '\%' << flush;
}


class ExonData{
public:
    t_exon start, end, five_to_three, three_to_five;

    ExonData(t_exon start_int, t_exon end_int, t_exon f_t_t, t_exon t_t_f){
        start = start_int;
        end = end_int;
        five_to_three = f_t_t;
        three_to_five = t_t_f;
    }
};


/** Reference object that each exon (GeneContainer) points to **/

class GeneStats{

//#define DEBUGGENE "ADAM15-ISOF2"
#define DEBUGGENE "PLA2R1-ISOF1"

#define debuggene(prefix, name,exon,start,stop,ft3,tt5) \
 if (name.startsWith(DEBUGGENE))\
    fprintf(stderr, "%s  %s Ex%02d [%9d,%9d]  (5>3) %6d, (3>5) %6d\n", prefix, name.toUtf8().data(), exon, start, stop, ft3, tt5); \
 else (void)0


public:
    QString gene;

    QMap<t_exon,ExonData *> exon_positions;
    // [#Exon] --> [
      //  pos1, pos2,
      //  cumulative coding length from 5->3,
      //  cumulative from 3->5]
    //

    GeneStats(QString gene){ this->gene = gene; }

    void insertExon(t_exon exon_num, t_exon pos1, t_exon pos2)
    {

        bool forward = true;
        if ((exon_num > 1 && !exonExists(exon_num-1)) || (exon_num <= 1 && exon_positions.size()>1)) forward = false;

        if (exonExists(exon_num)){
//            cerr << "Exon already exists!" << gene.toUtf8().data() << "--" << exon_num << endl;
        } else {
            t_exon diff= pos2 - pos1;
            t_exon rollover = diff;

//            cerr << this->gene.toUtf8().data() << endl;
//            cerr << exon_num << endl;

           //Carry on count for 5->3'
            t_exon previous_exon_number = exon_num + (forward?-1:+1);

            if (exonExists(exon_num-1)){
                t_exon previous_diff = exon_positions.value(previous_exon_number)->five_to_three;
                rollover += previous_diff;
            }

            //Update the cumulative 3->5 field for each exon below this
            t_exon prev_index = previous_exon_number;

            while (exonExists(prev_index)){
              exon_positions.value(prev_index)->three_to_five += diff;
              forward?prev_index --:prev_index++;
            }

            ExonData *data = new ExonData(pos1, pos2, rollover, diff);
//            data->start = pos1;
//            data->end = pos2;
//            data->five_to_three = rollover;
//            data->three_to_five = diff;

            debuggene("-->Inserted:",gene, exon_num, pos1, pos2, rollover, diff);
            exon_positions[exon_num] = data; //stores pointer

            //DEBUG: Works for forward encoded genes
            t_exon start_index = 0;
            while (!exonExists(++start_index) && start_index <= 50 );

            do {
                ExonData *prev = exon_positions.value(start_index);
                debuggene("exExists:",gene, start_index, prev->start, prev->end, prev->five_to_three, prev->three_to_five);
            } while (exonExists(++start_index));
        }
    }

    bool exonExists(t_exon exon_number){
        return exon_positions.contains(exon_number);
    }
};




/** Stores data about each exon **/
class GeneContainer{
public:
    QString name;
    QString chr;
    int pos1, pos2;
    int frame;
    bool forward;

    GeneContainer(QString &name, QString chr, int pos1, int pos2, int frame, bool forward){
        this->name = name;
        this->chr = chr;
        this->pos1 = pos1; this->pos2 = pos2;
        this->frame = frame;
        this->forward = forward;
    }
};



class GeneMap{
public:
    QMap<QString,QMap<QString,GeneContainer*> > Gene_Map;
    QMap<QString,QMap<QString,GeneStats*> > Gene_Stats;
    // chr1 --> Gene --> Container

    GeneMap(QString gmp_file){ populateGeneMap(gmp_file);}

private:
    void populateGeneMap(QString &gmp_file){

        uint numlines = countlines(gmp_file);
        uint cnt = 0;
        QFile inputFile(gmp_file);

        if (inputFile.open(QIODevice::ReadOnly))
        {
            QTextStream in(&inputFile);
            while (!in.atEnd())
            {
                QStringList tokens = in.readLine().split('\t');
                int len = tokens.length();

                if (len<6){
                    cerr << "genemap must contain reading frames and strand direction!" << endl;
                    exit(-1);
                }

                QString name = tokens[3].trimmed();
                QString chr = tokens[0].trimmed();

                t_exon pos1 = tokens[1].toInt();
                t_exon pos2 = tokens[2].toInt();

                GeneContainer *g = new GeneContainer(
                        name,                        // Name
                        chr,                        // Chr
                        pos1, pos2,                // Pos1, Pos2
                        tokens[len-1].toInt(),    // Frame
                        tokens[len-2][0]=='+'    // Forward
                );
                Gene_Map[chr][name] = g;


                // Determine the length of the coding sequence for this gene
                bool splice=false, exon=false;
                int Exon_number=-1;

                QStringList name_exon_splice = name.split('|');
                if (name_exon_splice.length()>1)
                {
                    QStringList exon_splice =  name_exon_splice[1].split('_');
                    Exon_number = exon_splice[0].split("Exon")[1].toInt();

                    if (exon_splice.length()>1)
                    {
                        if (exon_splice[0].trimmed().startsWith("Exon")) exon=true;

                        for(int si=1; si< exon_splice.size(); si++){
                            if (exon_splice[si].trimmed().startsWith("Splice")) {
                                splice=true; break;
                            }
                        }
                    }
                    else exon=true;
                } else exon = true; // single gene name

//                cerr << "gene=" << name.toUtf8().data() << ", exon=" << exon << ", splice=" << splice << endl;


                if (exon && !splice){
                    // Update GeneCDS positions
                    QString nameref = name_exon_splice[0].trimmed();

                    if (Gene_Stats[chr].contains(nameref)){
                        Gene_Stats[chr][nameref]->insertExon(Exon_number, pos1, pos2);
                    } else {
                        (Gene_Stats[chr][nameref] = new GeneStats(nameref))->insertExon(Exon_number, pos1, pos2);
                    }
                }
                if (++cnt%15321==0) progress("Map", 0, 100*cnt/numlines);
            }

            //THAR BE SPACE BUGS HERE YAR!
//            GeneStats *ree = Gene_Stats["chr1"]["HRNR"];
//            cerr << '\n' << ree->gene.toUtf8().data() << ", "
//                 << ree->exon_positions[2][0] << "-" << ree->exon_positions[2][1]
//                 << "==" << ree->exon_positions[2][2] << endl;
//            exit(-1);
            progress("Map", 0, 100);
            inputFile.close();
        }
    }

};


#endif // GENEMAP_H