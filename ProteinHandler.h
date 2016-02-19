#ifndef PROTEINHANDLER_H
#define PROTEINHANDLER_H

#define SYNONYMOUS "SYNONYMOUS"
#define MISSENSE   "MISSENSE"
#define NONSENSE   "NONSENSE"

#include "FASTAhandler.h"

#include <QMap>
#include <QStringList>

#include <iostream>

using namespace std;


/*
Map format:
Alu/A   AAC,TAC,ACT
Str/S   GER,GRE,VER
*/


class ProteinHandler{

    QMap<QString, QString> symbol_map;  //[ALA] -> A  keep it as QString for simplicity with other functions
    QMap<QChar,QChar> rcMap; //reverse complement map;

    void populateRCMap(){
        rcMap['A'] = 'T';rcMap['T'] = 'A';rcMap['a'] = 't';rcMap['t'] = 'a';
        rcMap['C'] = 'G';rcMap['G'] = 'C';rcMap['c'] = 'g';rcMap['g'] = 'c';
    }

    QStringList codon_list;
    void populateCodonMap(QString &dna_codon_file){
        QStringList codon_list;
        QStringList protein_list;
        QStringList symbol_list;

        QFile inputFile(dna_codon_file);

        if (inputFile.open(QIODevice::ReadOnly))
        {
            QTextStream in(&inputFile);
            while (!in.atEnd())
            {
                QStringList tokes = in.readLine().split('\t');
                QStringList protein_symb = tokes[0].trimmed().split('/');

                QString protein = protein_symb[0].trimmed();
                QString symb = protein_symb[1].trimmed();

                QStringList codons = tokes[1].trimmed().split(',');

                for (int c=0; c < codons.size(); c++){
                    QString codon = codons[c].trimmed();
                    if (!codon_list.contains(codon)){
                        codon_list.append(codon);
                    }  else {
                        cerr << "Error: Duplicate codon in DNA table: " << codon.toUtf8().data() << endl;
                        exit(-1);
                    }
                    protein_map[codon]=protein;
                }

                if (!protein_list.contains(protein)){
                    protein_list.append(protein);
                } else {
                    cerr << "Error: Duplicate protein name in DNA table" << protein.toUtf8().data() << endl;
                    exit(-1);
                }

                if (!symbol_list.contains(symb)){
                    symbol_list.append(symb);
                } else {
                    cerr << "Error: Duplicate protein symbol in DNA table" << symb.toUtf8().data() << endl;
                    exit(-1);
                }
                symbol_map[protein]=symb;
            }
            inputFile.close();
        }
    }


public:
    QMap<QString, QString> protein_map;

    ProteinHandler(QString dna_codon_file){
        populateCodonMap(dna_codon_file);
        populateRCMap();
    }


    QString codonToprotein(QString codon, bool letter=false){
        codon = codon.trimmed();

        QString codon_upper = codon.toUpper();
        QChar defval='@';

        QString protein = protein_map.value(codon_upper,defval);
        if (letter) protein = symbol_map.value(protein,defval);

        if (protein == defval){
            cerr << "Could not find key1= "
                 << codon_upper.toUtf8().data() << endl;
            exit(-1);
        }
        return protein;
    }


    QString getProteinSymbol(QString protein){
        QChar defval='@';
        QString symbol = symbol_map.value(protein,defval);
        if (symbol!=defval) return symbol;

        cerr << "Could not find key for:" << protein.toUtf8().data() << endl;
        exit(-1);
    }

    QString codons2Proteins(QString REF_codon, QString ALT_codon, bool letter=false){
        REF_codon = REF_codon.trimmed();
        ALT_codon = ALT_codon.trimmed();

        QString REF_upper = REF_codon.toUpper(), ALT_upper = ALT_codon.toUpper();
        QChar defval='@';
//        bool upper = true; //Assume all are upper

//        if (REF_codon != REF_upper) upper=false; //except when proven not
//        if (ALT_codon != ALT_upper) upper=false;

        QString ref_protein = protein_map.value(REF_upper,defval);
        QString alt_protein = protein_map.value(ALT_upper,defval);
        if (letter) {
            ref_protein = symbol_map.value(ref_protein,defval);
            alt_protein = symbol_map.value(alt_protein,defval);
        }

        if ((ref_protein == defval) || (alt_protein == defval)){
            cerr << "Could not find key2= "
                 << REF_codon.toUtf8().data() << ','
                 << ALT_codon.toUtf8().data() << endl;
            exit(-1);
        }
//        return ref_protein+','+alt_protein+','+(upper?'Y':'N');
        return ref_protein+','+alt_protein;
    }

    QString reverseComplement(QString codon){
        //three letters -- direct assignment is quicker.
        return QString()\
                .append(rcMap[codon[2]])\
                .append(rcMap[codon[1]])\
                .append(rcMap[codon[0]]);
    }

    const char * proteins2mutation(QString ref_protein, QString alt_protein){
        if (ref_protein==alt_protein) return SYNONYMOUS;
        if (alt_protein=="STOP") return NONSENSE;
        return MISSENSE;
    }


};


#endif // PROTEINHANDLER_H
