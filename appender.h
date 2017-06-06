#ifndef APPENDER_H
#define APPENDER_H

#include "GeneMap.h"


/* We include two type of headers:
 * For each gene/exon in a genelist, print a corresponding:
 *     + RAcodon  ([AAC ACA],[CCD DCD], )
 *     + RAcodchange ( g.76A>C, g.78C>D )   // genomic sequence change from codon list
 *     + RAprotein  ([A L],[R T], )
 *     + RAProchange ( p.A14L, p.R12T
 *     + mutation (MISS,NONS,)
where RA prints both the Reference and Alternate codon and protein.
Mutation is simply how the reference protein maps to the alt, so this data does not need to be paired.

 Thus we need 3 headers:
 */
#define C_ID "COL"    //Codon List
#define P_ID "PRL"    //Protein List
#define M_ID "MUL"    //Mutation List
#define D_ID  "DIL"    //Direction List

#define CH_ID "CCH"   //Codon change list   --  //Forward coordindates only
#define PH_ID "PCH"   //Protein change list --  //forwards coordinates only


#define HFORMAT "##FORMAT="
#define FILLER ",Number=.,Type=String,Description=\""
#define CAP "\">"

#define HEADER_CLIST HFORMAT "<ID=" C_ID
#define HEADER_PLIST HFORMAT "<ID=" P_ID
#define HEADER_MLIST HFORMAT "<ID=" M_ID
#define HEADER_DLIST HFORMAT "<ID=" D_ID

#define HEADER_COCHLIST HFORMAT "<ID=" CH_ID
#define HEADER_PRCHLIST HFORMAT "<ID=" PH_ID

#define HEADER_CLIST_FULL HEADER_CLIST FILLER "Codon list corresponding to a genelist" CAP
#define HEADER_PLIST_FULL HEADER_PLIST FILLER "Protein list corresponding to a genelist" CAP
#define HEADER_MLIST_FULL HEADER_MLIST FILLER "Mutation list corresponding to a genelist" CAP
#define HEADER_DLIST_FULL HEADER_DLIST FILLER "Direction list corresponding to a genelist" CAP

#define HEADER_COCHLIST_FULL HEADER_COCHLIST FILLER "Codon transition list corresponding to a genelist" CAP
#define HEADER_PRCHLIST_FULL HEADER_PRCHLIST FILLER "Protein transition list corresponding to a genelist" CAP

//Ala --> Val  these are harmless mutations
//Val --> Ala
//Val --> Gly
//Gly --> Val




class Appender{
private:
    //Writing methods
    void writeout(const char *message, bool *newline, QFile *&file){
        file->write(message);
        if (*newline) file->write("\n");
        file->flush();
    }

    void rejout(QString message, bool newline=true){ // assume string not in mem yet
        writeout(message.toUtf8().data(), &newline, rejects);
    }

    void fileout(const char * message, bool newline=true){
        writeout(message, &newline, outputs);
    } 


public:
    //indexes
    int INFO_INDEX, FORMAT_INDEX, INDIV_START_INDEX;  // updated by handleHeaders()
    int GENES_INDEX, TYP_INDEX; // updated every chromosome (should be the same though)

    //Static shared over all instances:/// HHHOWW?
    QMap<QString,QMap<QString,GeneContainer*> > genemap;        //generate once
    QMap<QString,QMap<QString,GeneStats*> > genestats;        //generate once

    ProteinHandler *ph;                                         //generate once
    FASTAHandler *fh;

    QFile *outputs;
    QFile *rejects;

    int tab_count;

    void handleHeaders(QString &line, QTextStream &in, uint &countline);

    bool isaSNV( QString &info_text);


    QString getRefCodon(bool &forward, int &ref_index, quint64 &bpos);
    QString getAltCodon(bool &forward, int &ref_index, bool &isSNV, bool &isDel,
                                  QString &ref_codon, QString &VALT,
                                  int &vrlen, int &valen, quint64 &bpos);

    QString c_nomenclature(int &position, QString &ref, QString &alt, bool &isSNV, bool &isDel);
    QString p_nomenclature(quint64 &bpos, int &position, QString &ref_p, QString &alt_protein, int vrlen, int valen, bool &isSNV, bool &isDel, quint64 &next_codon_bpos);

    QStringList getLists(QString &chr, quint64 &bpos, QStringList &genelist,
                              QString &VALT, QString &VREF, int vrlen, int valen, bool &isSNV, bool &isDel, QString &rejects_per_line);

    QString antonorakis(ExonData *data, bool &direction, quint64 &bpos, QString &VREF, QString &VALT,
                                  QString &ref_codon, QString &alt_codon,
                                  bool &isSNV, bool &isDel , quint64 &next_codon_bpos, bool printProtein);

    ~Appender(){
        outputs->close();
        rejects->close();
        delete fh;
    }

    Appender(QString &vcf_file, QString &FASTA_folder, QString &G_id,
             QString &outfold, QString &rejfold,
             int tabcount,
             GeneMap &map, ProteinHandler &prh)
    {
        QString prefix = ((vcf_file.contains('/'))?vcf_file.split('/').last():vcf_file);
        prefix = (prefix.endsWith(".vcf")?prefix.split(".vcf")[0]:prefix);

        const char * prefs = prefix.toUtf8().data();

        cerr << prefs << ":  " << flush;

        (outputs = new QFile(outfold+'/'+prefix+".func.vcf"))->open(QIODevice::WriteOnly);
        (rejects = new QFile(rejfold+'/'+prefix+".func.rejects"))->open(QIODevice::WriteOnly);

        //non-s
        fh = new FASTAHandler(FASTA_folder);

        // Would be static if it was possible to assign references to a class that may/may not be instantiated (virtual?)
        // Could create a static holder parent -- elegant but wasteful. Just pass pointers each time -- easier.
        ph = &prh;
        genemap = map.Gene_Map;
        genestats = map.Gene_Stats;

        QFile inputFile(vcf_file);

        uint numlines = countlines(vcf_file);
        uint countline = 0;

        //Tmp variables
        QString current_chr="";

        if (inputFile.open(QIODevice::ReadOnly))
        {
            QTextStream in(&inputFile);
            QString line="##";
            handleHeaders(line, in, countline);

            while (!in.atEnd())
            {
//                cerr << line.toUtf8().data() << endl;


                //FILE number (i.e. first, second, third ,etc) corresponds to num tabs -- do a macro for this
                // also rejects should appear in same line possibly with [sympton] -- DONE

                if(++countline%12==0) progress(100*countline/numlines);

                QString line = in.readLine();
                if (line.trimmed().length()<1) continue;
                QStringList tokes = line.split('\t');

                QString chr = tokes[0].trimmed();
                quint64 pos = tokes[1].toULongLong();

                QStringList FORMAT = tokes[FORMAT_INDEX].split(':');
                QStringList IFORMAT = tokes[INDIV_START_INDEX].split(':');

//                cerr << "IFORMAT == " << tokes[INDIV_START_INDEX].toUtf8().data() << endl;

                // Chromosomes dont match, update FASTA
                // update indexes too
                if (current_chr!=chr){
                    fh->openFASTA(chr);
                    current_chr = chr;

                   GENES_INDEX = FORMAT.indexOf(G_id);
                }

                if (GENES_INDEX == -1){
                    cerr << "\n[Error] Could not find " << G_id.toUtf8().data();
                    cerr << " for line " << (countline+1);
                    cerr << " in file: " << vcf_file.toUtf8().data() << endl;

                    cerr << "GI=" << GENES_INDEX << endl;

                    exit(-1);
                }
                //Add id's to FORMAT
                FORMAT.insert(GENES_INDEX+1,PH_ID);
                FORMAT.insert(GENES_INDEX+1,CH_ID);
                FORMAT.insert(GENES_INDEX+1,M_ID);
                FORMAT.insert(GENES_INDEX+1,P_ID);
                FORMAT.insert(GENES_INDEX+1,C_ID);
                FORMAT.insert(GENES_INDEX+1,D_ID);

                QString newformat;
                for (int t=0; t< FORMAT.length(); t++) newformat.append(FORMAT[t]).append(':');
                newformat.chop(1); // remove last colon
                tokes[FORMAT_INDEX] = newformat;


                //Grab VCF reference, VCF alt, FASTA reference
                QString VREF = tokes[3].trimmed();
                QString VALT = tokes[4].trimmed();
                QString FREF = fh->getReference(pos);

                bool isSNV = isaSNV(tokes[INFO_INDEX]);
                bool isDel = (!isSNV && VALT.length()==1);

//                if (!isSNV){
//                    cerr << "\nNot a SNV:\n" << chr.toUtf8().data() << "  " << pos << endl;
//                    cerr << IFORMAT[TYP_INDEX].trimmed().toUtf8().data() << endl;
//                    exit(-1);
//                }



                QString rejects_per_line=" "; //collects bad genes/introns/references etc
                if (VREF[0]!=FREF[0]){
                    rejects_per_line.append("[VCF_REF="+VREF
                                            +" FASTA_REF="+FREF
                                            +" differ in first letter, converting to upper and ignoring],");
                }
                //directions, codons, proteins, mutations

//                cerr << "gettinglists:- " << countline << flush;
//                cerr << line.toUtf8().data() << flush;
//                cerr << "isSNV, isDel: " << isSNV << "," << isDel << flush;

                QStringList genelist = IFORMAT[GENES_INDEX].split(',');
//                cerr << chr.toUtf8().data() << " " << pos << " [" << flush;
//                for (int i=0; i< genelist.length(); i++) cerr << " " << genelist[i].toUtf8().data() << flush;
//                cerr << "] " << VALT.toUtf8().data() << "=" << VREF.toUtf8().data() << " (" << VREF.length() << "," << VALT.length() << ")" << isSNV << " " << isDel << " " << rejects_per_line.toUtf8().data() << endl;

                QStringList d_c_p_m_ch_ph = getLists(chr, pos, genelist, VALT, VREF, VREF.length(), VALT.length(), isSNV, isDel, rejects_per_line);

//                cerr << " --> gotlists" << endl;
//                cerr << '\n' << line.toUtf8().data() << endl;

                //Add lists to IFORMAT
                IFORMAT.insert(GENES_INDEX+1,d_c_p_m_ch_ph[5]);
                IFORMAT.insert(GENES_INDEX+1,d_c_p_m_ch_ph[4]);
                IFORMAT.insert(GENES_INDEX+1,d_c_p_m_ch_ph[3]);
                IFORMAT.insert(GENES_INDEX+1,d_c_p_m_ch_ph[2]);
                IFORMAT.insert(GENES_INDEX+1,d_c_p_m_ch_ph[1]);
                IFORMAT.insert(GENES_INDEX+1,d_c_p_m_ch_ph[0]);
                QString newiformat;
                for (int t=0; t< IFORMAT.length(); t++) newiformat.append(IFORMAT[t]).append(':');
                newiformat.chop(1);
                tokes[INDIV_START_INDEX] = newiformat;

                //print line
                int last = tokes.length()-1;
                for (int ty=0; ty <= last; ty++){
                    fileout(tokes[ty].trimmed().toUtf8().data(), false);
                    if (ty!=last) fileout("\t",false);
                }

                //Handle rejects
                rejects_per_line.chop(1);
                if(rejects_per_line.length()>1) rejout(chr+'\t'+QString::number(pos)+'\t'+rejects_per_line);

                fileout("\n",false);
            }
            progress(100);
            cerr << endl;
        }
    }
};

#endif // APPENDER_H
