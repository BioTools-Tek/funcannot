#ifndef FASTAHANDLER_H
#define FASTAHANDLER_H

#include <QFile>
#include <QTextStream>
#include <iostream>

#define sym_per_row 50
#define chr_per_row 51

class FASTAHandler{
public:
    QString folder;
    uint offset;
    QFile *FAST_FILE;
    QTextStream *FASTA;

    FASTAHandler(QString fast_folder){
        //Set null pointers + defaults
        FAST_FILE = 0; FASTA = 0;
        folder = fast_folder;
    }

    void openFASTA(QString chr){
        //Close old FASTA
        if (FAST_FILE!=0){
            if (FAST_FILE->isOpen()) FAST_FILE->close();
        }

        //Open new, update pointers
        FAST_FILE = new QFile(this->folder+'/'+chr+".fa");
        if (FAST_FILE->open(QIODevice::ReadOnly)) FASTA = new QTextStream(FAST_FILE);

        offset = FASTA->readLine().length()+1; // +1 for newline
    }


    QString getReference(quint64 bp, int length = 1, bool seekbackfirstchar=false, bool printRC = false )
    {
        quint64 bpos = bp - 1;
        quint64 rows = bpos/sym_per_row;
        quint64 base_row = rows*sym_per_row;
        uint cols = bpos - base_row;

        if (printRC)
            return QString("[row=").append(QString::number(rows+2)).append(" col=").append(QString::number(cols)).append(']');

        quint64 chars = offset + base_row + rows + cols;
        if (FAST_FILE->isOpen())
            FASTA->seek(chars);
        else {
            std::cerr << "Cannot open FASTA:" << FAST_FILE->fileName().toUtf8().data() << "  " << bp << std::endl;
            return "";
        }
        // the default FASTA->read(3) may contain newlines and introns, need to skip

        QString res="";

        bool first_char=true;
        // If the first character is lower case, move backwards until it finds an uppercase one and then resume search
        // but only if seekbackfirstchar is true.

        while (length > 0){
            QChar base = FASTA->read(1)[0];
            //Skip new lines
            if (base=='\n') continue;

            //Skip introns (only want coding)
            QChar base_upper = base.toUpper();
            if (base!=base_upper) {
                if (seekbackfirstchar && first_char){
                    chars -= 2; //move back one character (2 to counter each read which shift pos forward)
                    FASTA->seek(chars);
                }
                continue;
            }
            res.append(base);
            length --;

            first_char = false;
        }

        return res;
    }
};





#endif // FASTAHANDLER_H
