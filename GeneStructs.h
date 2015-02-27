#ifndef GENESTRUCTS_H
#define GENESTRUCTS_H

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


#ifdef DEBUG
    #define DEBUGGENE "ADAM15-ISOF2"
    #define DEBUGGENE "PLA2R1-ISOF1"

    #define debuggene(prefix, name,exon,start,stop, diff, ft3,tt5) \
 if (name.startsWith(DEBUGGENE))\
    fprintf(stderr, "%s %s Ex%02d [%9d-%9d =%4d] %5d(5>3) %5d(3>5)\n", prefix, name.toUtf8().data(), exon, start, stop, diff, ft3, tt5); \
 else (void)0

#endif

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
            const t_exon diff= pos2 - pos1;
            t_exon rollover = diff;

//            cerr << this->gene.toUtf8().data() << endl;
//            cerr << exon_num << endl;

           //Carry on count for 5->3'
            t_exon previous_exon_number = exon_num + (forward?-1:+1);

            if (exonExists(previous_exon_number)){
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

#ifdef DEBUG
            debuggene("-->Inserted:",gene, exon_num, pos1, pos2, diff, rollover, diff);
#endif
            exon_positions[exon_num] = data; //stores pointer

#ifdef DEBUG
            //DEBUG: Works for forward encoded genes
            t_exon start_index = 0;
            while (!exonExists(++start_index) && start_index <= 50 );

            do {
                ExonData *prev = exon_positions.value(start_index);
                debuggene("exExists:",gene, start_index, prev->start, prev->end, prev->end-prev->start,  prev->five_to_three, prev->three_to_five);
            } while (exonExists(++start_index));
#endif
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

#endif // GENESTRUCTS_H
