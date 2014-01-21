#ifndef GENEMAP_H
#define GENEMAP_H

#include "ProteinHandler.h"

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


/** Reference object that each exon (GeneContainer) points to **/
class GeneStats{
public:
    QString gene;

    QMap<int,QList<int> > exon_positions;
    // [#Exon] --> [pos1, pos2, cumulative coding length from 5->3, cumulative from 3->5]

    GeneStats(QString gene){ this->gene = gene; }

    void insertExon(int exon_num, int pos1, int pos2){
        if (exonExists(exon_num)){
            cerr << "Exon already exists!" << gene.toUtf8().data() << "--" << exon_num << endl;
        } else {
            int diff= pos2 - pos1;
            int rollover = diff;

            //Carry on count for 5->3'
            if (exonExists(exon_num-1)){
                int previous_diff = exon_positions[exon_num-1][2];
                rollover += previous_diff;
            }

            //Update the cumulative 3->5 field for each exon below this
            int prev_index = exon_num-1;
            while (exonExists(prev_index)){
                exon_positions[prev_index][3] += diff;
                prev_index --;
            }

            QList<int> data;
            data.append(pos1);
            data.append(pos2);
            data.append(rollover); //cumulative total, including the current exon
            data.append(diff);
            exon_positions[exon_num] = data;
        }
    }

    bool exonExists(int exon_number){
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

                int pos1 = tokens[1].toInt(), pos2 = tokens[2].toInt();

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
                            if (exon_splice[1].trimmed().startsWith("Splice")) splice=true;
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
