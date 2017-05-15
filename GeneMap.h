#ifndef GENEMAP_H
#define GENEMAP_H

//#define DEBUG

#include "GeneStructs.h"


// Used to parse headers in GeneMap
#define HEAD_CHROM "#CHROM"
#define HEAD_START "START"
#define HEAD_STOP  "STOP"
#define HEAD_GINF  "GENEINFO"
#define HEAD_SCORE "SCORE"          //Score we don't really care about, handled upstream
#define HEAD_DIRECT "DIRECT"
#define HEAD_NMINFO "RGNAME"        //This is ignored
#define HEAD_FRAMES "FRAMES"


class GeneMap{
public:
    QMap<QString,QMap<QString,GeneContainer*> > Gene_Map;
    QMap<QString,QMap<QString,GeneStats*> > Gene_Stats;
    // chr1 --> Gene --> Container

    GeneMap(QString gmp_file){ populateGeneMap(gmp_file);}

private:
    QMap <QString, int> column_map;

    void mapHeaders(QString line){
        if (!line.startsWith('#')){
            cerr << "A header line was not found in the gene map." << endl;
            exit(-1);
        }

        QStringList tokes = line.split('\t');

        for (int i=0; i< tokes.length(); i++){
            QString t = tokes[i].trimmed();
            column_map[t] = i;
        }

        //Test essentials
        if (column_map[HEAD_CHROM] == -1 || column_map[HEAD_START] == -1 || column_map[HEAD_STOP] == -1 || column_map[HEAD_GINF] == -1){
            cerr << "Either chrom, start, stop, or geneinfo are not present in genemap. Quitting." << endl;
            exit(-1);
        }

        //Test for reading frame and direction
        if (column_map[HEAD_FRAMES] == -1 || column_map[HEAD_DIRECT] == -1){
            cerr << "Genemap must contain reading frames and strand direction!" << endl;
            exit(-1);
        }
    }


    void populateGeneMap(QString &gmp_file){

        uint numlines = countlines(gmp_file);
        uint cnt = 0;
        QFile inputFile(gmp_file);

        if (inputFile.open(QIODevice::ReadOnly))
        {
            QTextStream in(&inputFile);

            //Process headers --> populate column map
            mapHeaders(in.readLine());

            cerr << "-Map:" << flush;

            while (!in.atEnd())
            {
                QStringList tokens = in.readLine().split('\t');

                QString chr  = tokens[column_map[HEAD_CHROM]].trimmed();
                t_exon pos1  = tokens[column_map[HEAD_START]].toInt(),
                       pos2  = tokens[column_map[HEAD_STOP]].toInt();

                QString name = tokens[column_map[HEAD_GINF]].trimmed();

                int frame    = tokens[column_map[HEAD_FRAMES]].toInt();
                bool forward = tokens[column_map[HEAD_DIRECT]][0] == '+';


                GeneContainer *g = new GeneContainer(
                        name,                        // Name
                        chr,                        // Chr
                        pos1, pos2,                // Pos1, Pos2
                        frame,                    // Frame
                        forward                  // Forward
                );
                Gene_Map[chr][name] = g;


                // Determine the length of the coding sequence for this gene
                bool splice=false, exon=false;
                int Exon_number=-1;

                QStringList name_exon_splice = name.split('|');
                if (name_exon_splice.length()>1)
                {
                    QStringList exon_splice =  name_exon_splice[1].split('_');

                    // Skip non-coding
                    if (!exon_splice[0].startsWith("Exon")){
                        continue;
                    }

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
                if (++cnt%15321==0) progress(100*cnt/numlines);
            }

            //THAR BE SPACE BUGS HERE YAR!
//            GeneStats *ree = Gene_Stats["chr1"]["HRNR"];
//            cerr << '\n' << ree->gene.toUtf8().data() << ", "
//                 << ree->exon_positions[2][0] << "-" << ree->exon_positions[2][1]
//                 << "==" << ree->exon_positions[2][2] << endl;
//            exit(-1);
            progress(100);
            cerr << endl;

            inputFile.close();
        }
    }

};


#endif // GENEMAP_H
