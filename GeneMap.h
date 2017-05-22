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


                QStringList gene_bar_detail = name.split('|');
                // [Name]                  '|'   '_'
                // Gene                   - 1  -  0
                // Gene|Exon1             - 2  -  0  * only scenario we have to work on
                // Gene|Exon1_SpliceA     - 2  -  2
                // Gene|Exon1_5'UTR       - 2  -  2
                // Gene|Intron1           - 2  -  0
                // Gene|Promoter          - 2  -  0
                // Gene1-Gene2|Intergenic - 2  -  0

                // However, this is FUNCTIONAL annotation -- meaning we only care about coding.

                if (gene_bar_detail.length()>1)
                {
                    QString gene_part   = gene_bar_detail[0]; // Gene
                    QString detail_part = gene_bar_detail[1]; // Exon1_5'UTR

                    /*
                    // Intergenic
                    if (detail_part.startsWith("Intergenic")){
                        //mapIntergenic();
                        continue;
                    }

                    // Promoter
                    if (detail_part.startsWith("Promoter")){
                        //mapPromoter();
                        continue;
                    }

                    // Intron
                    if (detail_part.startsWith("Intron")){
                        continue;
                    }*/

                    if (detail_part.startsWith("Exon")){

                        // Pure coding region -- not splice, not UTR, boom we're in business.
                        if (!detail_part.contains('_')){
                            int exon_number = detail_part.split("Exon")[1].toInt();

                            // Update GeneCDS positions
                            QString nameref = gene_part.trimmed();

                            if (Gene_Stats[chr].contains(nameref)){
                                Gene_Stats[chr][nameref]->insertExon(exon_number, pos1, pos2);
                            } else {
                                (Gene_Stats[chr][nameref] = new GeneStats(nameref))->insertExon(exon_number, pos1, pos2);
                            }
                        }
                    }
                }
                if (++cnt%15321==0) progress(100*cnt/numlines);
            }

            progress(100);
            cerr << endl;

            inputFile.close();
        }
    }

};


#endif // GENEMAP_H
