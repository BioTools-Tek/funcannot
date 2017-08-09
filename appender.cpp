#include "appender.h"

#define BSta '['
#define BEnd ']'
#define badId "*,"


#define skipbad {\
    protein_change_list.append(badId);\
    codon_change_list.append(badId);\
    codonlist.append(badId);\
    proteinlist.append(badId);\
    mutationlist.append(badId);\
    directionlist.append(badId);\
    continue;\
}


// Ref_index is in relation to direction. 0 - first part of codon, 1 - second part of codon, 2 - third part of codon in both cases
// reverse positions require a +1 offset it seems (ref_index, getRefCodon) ... not sure why.

QString Appender::getRefCodon(bool &forward, int &ref_index, quint64 &bpos){
//    cerr << "getRefCodon: " << "forw=" << forward << " ref=" << ref_index << " bpos=" << bpos << endl;
    if(forward) return fh->getReference(bpos-ref_index,3);
    return fh->getReference(bpos - (3 - ref_index) + 1,3); // No reverse complimenting done here. getAltCodon relies on unaltered ref.
}

QString Appender::getAltCodon(bool &forward, int &ref_index, bool &isSNV, bool &isDel,
                              QString &ref_codon, QString &VALT,
                              int &vrlen, int &valen, quint64 &bpos){

    //Possible scenarios:                                                              Key: S = isSNV, D = isDel, F = Forward, rN = ref_index equals N
    //                                                                                      lowercase negates booleans (SDF)
// 1a. SNV change on the first index(+)/ last index(-)
    // ref_index:       0               2
    // VALT:            C               C
    //                A|ATG|C         A|ATG|C
    // direction:  >=========>      <========<
    // res:            C+TG         reverseCompl(C+TG)
    // params:       SdFr0           Sdfr2

// 1b. SNV change on middle index (+)/(-)
    // ref_index:        1               1
    // VALT:             C               C
    //                A|ATG|C         A|ATG|C
    // direction:  >=========>      <========<
    // res:             A+C+G       reverseCompl(A+C+G)
    // params:        SdFr1            Sdfr1

// 1c. SNV change on last index (+) / first index (-)
    // ref_index:        2               0
    // VALT:              C               C
    //                A|ATG|C         A|ATG|C
    // direction:  >=========>      <========<
    // res:             AT+C       reverseCompl(AT+C)
    // params:        SdFr2           Sdfr0


//    cerr << "getAltCodon: forw=" << forward << " ref=" << ref_index << " SNV=" << isSNV << " DEL=" << isDel
//         << " ref_codon=" << ref_codon.toUtf8().data() << " VALT=" << VALT.toUtf8().data() << " vrlen=" << vrlen << " valen=" << valen
//         << " bpos=" << bpos << endl;

    if (isSNV){
        if (forward){
            if (ref_index==0) return VALT+ref_codon.mid(1,2);                                           //1a +
            if (ref_index==1) return ref_codon.at(0)+VALT+ref_codon.at(2);                              //1b +
            if (ref_index==2) return ref_codon.left(2)+VALT.at(0);                                      //1c +
        }
        else {
            if (ref_index==0) return ph->reverseComplement(ref_codon.left(2)+VALT.at(0));               //1c -
            if (ref_index==1) return ph->reverseComplement(ref_codon.at(0)+VALT+ref_codon.at(2));       //1b -
            if (ref_index==2) return ph->reverseComplement(VALT+ref_codon.mid(1,2));                    //1a -
        }
    }

    //Deletions in VCF file are always given in the 5'->3' direction
// 2a. Deletion on first index(+) / last index(-)
    // ref_index:       0                      2
    // VREF->VALT:     AATGC -> A          AATGC -> A
    //               T|AAT|GCA|T         T|AAT|GCA|T
    // direction:   >===========>       <===========<
    // res:            AAT                reverseCompl(AAT)
    // params:       sDFr0                  sDfr2

// 2b. Deletion on second index (+)/(-)
    // ref_index:       1                      1
    // VREF->VALT:      ATGC -> A           ATGC -> A
    //               T|AAT|GCA|T         T|AAT|GCA|T
    // direction:   >===========>       <===========<
    // res:            A+AT               reverseCompl(A+AT)
    // params:        sDFr1                 sDfr1

// 2c. Deletion on last index (+)/ first index (-)
    // ref_index:       2                      0
    // VREF->VALT:       TGC -> A            TGC -> A
    //               T|AAT|GCA|T         T|AAT|GCA|T
    // direction:   >===========>       <===========<
    // res:            AA+T               reverseCompl(AA+T)
    // params:        sDFr2                 sDfr0

    if (isDel){
        if (forward){
            if (ref_index==0) return (VALT + fh->getReference(bpos+vrlen,3)).left(3);                                           //2a +
            if (ref_index==1) return (ref_codon.at(0) + VALT + fh->getReference(bpos+vrlen,3)).left(3);                         //2b +
            if (ref_index==2) return ref_codon.left(2)+VALT.at(0);                                                              //2c +
        }
        else {
            if (ref_index==0) return ph->reverseComplement(ref_codon.left(2)+VALT.at(0));                                       //2c -
            if (ref_index==1) return ph->reverseComplement((ref_codon.at(0) + VALT + fh->getReference(bpos+vrlen,3)).left(3));  //2b -
            if (ref_index==2) return ph->reverseComplement((VALT + fh->getReference(bpos+vrlen,3)).left(3));                    //2a -
        }
    }

    //Insertions in VCF file are also given in in 5'->3' direction
// 3a. Insertion on the first index(+) / last index(-)
    // ref_index:         0                       2
    // VREF->VALT      A->ATCTC                A->ATCTC
    //              TT|ATG|TCC              TT|ATG|TCC
    // direction:  >==========>            <==========<
    // res:            ATC+                  reverseComp(ATC+)
    // params:       sdFr0                   sdfr2

// 3b. Insertion on the second index (+)/(-)
    // ref_index:         1                       1
    // VREF->VALT       T->TCTC                 T->TCTC
    //              TT|ATG|TCC              TT|ATG|TCC
    // direction:  >==========>            <==========<
    // res:            A+TT                  reverseComp(A+TT)
    // params:       sdFr1                   sdfr1

// 3c. Insertion on the last index (+)/ first index(-)
    // ref_index:         2                       0
    // VREF->VALT        G->TCTC                 G->TCTC
    //              TT|ATG|TCC              TT|ATG|TCC
    // direction:  >==========>            <==========<
    // res:            AT+T                  reverseComp(AT+T)
    // params:       sdFr2                   sdfr0

    if (forward){
        if (ref_index==0) return (VALT + fh->getReference(bpos+valen,2)).left(3);                                          // 3a +
        if (ref_index==1) return (ref_codon.at(0) + VALT + fh->getReference(bpos+1,2)).left(3);                            // 3b +
        if (ref_index==2) return ref_codon.left(2)+VALT.at(0);                                                             // 3c +/ 2c+/ 1c+
    }
    else {
        if (ref_index==0) return ph->reverseComplement(ref_codon.left(2)+VALT.at(0));                                       // 3c -/2c-/ 1c-
        if (ref_index==1) return ph->reverseComplement((ref_codon.at(0) + VALT + fh->getReference(bpos+1,2)).left(3));      // 3b -
        if (ref_index==2) return ph->reverseComplement((VALT + fh->getReference(bpos+valen,2)).left(3));                    // 3a -
    }

    cerr << "something went horrificly wrong" << endl;
    exit(-1);
}


//Constructer methods in header, as are indexes

QStringList Appender::getLists(QString &chr, quint64 &bpos, QStringList &genelist,
                          QString &VALT, QString &VREF, int vrlen, int valen,
                               bool &isSNV, bool &isDel,
                               QString &rejects_per_line){

    QString codonlist,         proteinlist,         mutationlist,        directionlist;
    QString codon_change_list, protein_change_list;

//    bpos -= 1;

    for (int g=0; g < genelist.length(); g++){
        QString gene = genelist[g].trimmed();
        QString just_gene= gene;
        t_exon exon_number=-1;

        int regulatory_val = Exon;  // Exon
        int regulatory_type = None; // Coding

        //Intron check
        if (gene.contains('|')){

            QStringList gene_extra = gene.split('|');

            just_gene = gene_extra[0];
            QString detail = gene_extra[1];

            // Skip non-exons
            if (detail.startsWith("Exon")){
                // Exons -- coding!
                QString num_extra = detail.split("Exon")[1];
                if(num_extra.contains("_")){
                    if (num_extra.contains("Splice")){
                        // Splice
                        QStringList deets = num_extra.split("_");

                        // Donor gives, Acceptor recieves
                        bool upstream = deets[1].split("Splice")[1] == "D";
                        num_extra = deets[0];

                        regulatory_val  = Splice_UTR;
                        regulatory_type = upstream?Upstream:Downstream;
                    }
                    else if (num_extra.contains("UTR")) {
                        rejects_per_line.append(BSta + gene + BEnd + ',');
                        skipbad;
                    }
                }
                exon_number = num_extra.toInt();
            }
            else {
                if (detail.startsWith("Promoter")){
                    bool upstream = detail.endsWith("upstream");

                    regulatory_val  = Promoter;
                    regulatory_type = upstream?Upstream:Downstream;
                }
                else {
                    rejects_per_line.append(BSta + gene + BEnd + ',');
                    skipbad;
                }
            }
        }

        if ( (!genemap.contains(chr)) || (!genemap[chr].contains(gene)) ){
            rejects_per_line.append(BSta +chr+'/'+gene+" not in genemap!"+ BEnd+',');
            skipbad;
        }

        GeneContainer *gc = genemap[chr][gene];
        GeneStats *gs = 0;

        if (genestats[chr].contains(just_gene)){
            gs = genestats[chr][just_gene];
        }
        else {
            //e.g. chr1 LOC100288069 -- complete UTR gene, promoter is therefore missed
            continue;
        }


        // If we assume that exonFrames refers to previous exon in terms of position, and not direction:

        // Discrepenacy between cDNA and Predicted Protein track.. two scenarios:
        // 1. cDNA is correct  // {0:sta, 1:mid, 2:end}/rev{0:end, 1: mid, 2:sta}
        // 2. Predicted Protein track is correct.
        // Assuming the cDNA, since more verifiable.
        int ref_index = (gc->forward)?\
((bpos - gc->pos1 + gc->frame - 1)%3):\
((gc->pos2 - bpos + gc->frame)%3);

        // 2. Predicted Protein track is correct.
//        int ref_index = (bpos - ( gc->pos1 + (3- gc->frame)) + 2) % 3;// {0:sta, 1:mid, 2:end}/rev{0:end, 1: mid, 2:sta}

        QString ref_codon = getRefCodon(gc->forward, ref_index, bpos);
        QString alt_codon = getAltCodon(gc->forward, ref_index, isSNV, isDel, ref_codon, VALT, vrlen, valen, bpos);

        quint64 next_ref_codon_index = bpos+(3-ref_index); // where the NEXT codon is on the reference

        //If reverse, then must complement:
        if (!gc->forward){
            directionlist.append("-,");
            ref_codon = ph->reverseComplement(ref_codon);  // -- This is needed here (can be reversed only AFTER getAltCodon has used it)
            VREF = ph->reverseComplement(VREF);
            VALT = ph->reverseComplement(VALT);
//            alt_codon = ph->reverseComplement(alt_codon);  // -- This is now performed in getAltCodon
        }
        else directionlist.append("+,");

        //This method counts from the start of the current exon
//        QStringList cnom_pnom_proteins = antonorakis(gs->exon_positions[Exon_number], gc->forward, bpos, VREF, VALT, ref_codon, alt_codon, isSNV, isDel, next_ref_codon_index ).split("||");

        //New method counts from start of coding sequence (first/last exon)

        ExonData * ex;

        if (regulatory_val == Promoter)
        {
            QList<t_exon> exon_keys = gs->exon_positions.keys();
            t_exon smallest = 9999999999999999999, largest = 0;
            t_exon smallest_value = 9999999999999999, largest_value = 0;

            for (int ek=0; ek < exon_keys.length(); ek ++){
                t_exon tmp = exon_keys.at(ek);

                t_exon start = gs->exon_positions.value(tmp)->start,
                       finit = gs->exon_positions.value(tmp)->end;

                t_exon small = (start < finit)?start:finit;
                t_exon large = (start > finit)?start:finit;


                if (small < smallest_value){
                    smallest = tmp;
                    smallest_value = small;
                }

                if (large > largest_value){
                    largest = tmp;
                    largest_value = large;
                }
            }

            ex = gs->exon_positions.value(gc->forward?smallest:largest);
        }

        /*else if (regulatory_val == Splice_UTR){
            ex = gs->exon_positions.value(exon_number);
        }*/

        else {
            // Regular Exon?
            ex = gs->exon_positions.value(exon_number);
        }

        QStringList cnom_pnom_proteins = antonorakis(
                    ex,
                    gc->forward,
                    bpos,
                    VREF, VALT,
                    ref_codon, alt_codon,
                    isSNV, isDel,
                    next_ref_codon_index, regulatory_type).split("||");


        if (regulatory_val == Exon){
            codon_change_list.append(cnom_pnom_proteins[0]+',');
            codonlist.append(BSta+ref_codon+' '+alt_codon+BEnd+',');

            QStringList proteins = cnom_pnom_proteins[2].split(' ');
            QString &ref_protein = proteins[0], &alt_protein = proteins[1];

            protein_change_list.append(cnom_pnom_proteins[1]+',');
            proteinlist.append(BSta+cnom_pnom_proteins[2]+BEnd+',');

            //Check Mutation
            QString mutation(ph->proteins2mutation(ref_protein, alt_protein));
            mutationlist.append(mutation+',');
        }
        else {
            codon_change_list.append(cnom_pnom_proteins[0]+',');
            codonlist.append("*,");

            protein_change_list.append("*,");
            proteinlist.append("*,");
            mutationlist.append("*,");
        }
    }

    QStringList res;
     // remove last comma
    codonlist.chop(1); proteinlist.chop(1); mutationlist.chop(1);
    directionlist.chop(1);codon_change_list.chop(1); protein_change_list.chop(1);

    res.append(directionlist); res.append(codonlist); res.append(proteinlist); res.append(mutationlist);
    res.append(codon_change_list);res.append(protein_change_list);
    return res;
}

QString Appender::c_nomenclature_regulatory(int &coding_pos, int &offset,
                                            QString &ref, QString &alt,
                                            bool &isSNV, bool &isDel)
{
    QString position = ((offset > 0)?"+":"")+QString::number(offset); // '-' in number

    if (coding_pos != 0){
        position.prepend(QString::number(coding_pos));
    }

    if (isSNV) return "c."+position+ref[0].toLatin1()+'>'+alt[0].toLatin1();
    if (isDel) return "c."+position+'_'+QString::number(coding_pos+ref.length())+"del"+ref;

    //insertion
    return "c."+position+'_'+QString::number(coding_pos+1)+"ins"+alt;
}


QString Appender::c_nomenclature(int &pos, QString &ref, QString &alt,
                              bool &isSNV, bool &isDel) {

    QString position = QString::number(pos);

    if (isSNV) return "c."+position+ref[0].toLatin1()+'>'+alt[0].toLatin1();
    if (isDel) return "c."+position+'_'+QString::number(pos+ref.length())+"del"+ref;

    //insertion
    return "c."+position+'_'+QString::number(pos+1)+"ins"+alt;
}


QString Appender::p_nomenclature(quint64 &bpos, int &position, QString &ref_protein, QString &alt_protein,
                              int vrlen, int valen,
                              bool &isSNV, bool &isDel,
                              quint64 &next_codon_bpos){

    QString ref_p = ph->getProteinSymbol(ref_protein),
            alt_p = ph->getProteinSymbol(alt_protein);

    if (isSNV) return "p."+ref_p+QString::number(position)+alt_p;
    if (isDel){
        int num_proteins_del = ((vrlen-1)/3)+1;
        if (num_proteins_del == 1) return "p."+ref_p+QString::number(position)+"del";
        return "p."+ref_p+QString::number(position)+'_'+alt_p+QString::number(position+num_proteins_del)+"del";
    }
    //else: insertion
    //<POS>A_<POS>BinsC
    //<POS>(ref_P)_<POS>(ref_P+1)ins<P up to and including alt_P>
    //Next_codon FASTA is only retrieved as and when needed -- i.e. here
    int num_proteins_ins = ((valen-1)/3)+1;

    QString intermediate_proteins="";

    quint64 next_three = next_codon_bpos-3;
//    cerr << "Pos,Next, valen pos: " << bpos << "," << next_three << ',' << valen << endl;
    while ( next_three < bpos+valen){
        intermediate_proteins.append( ph->codonToprotein(fh->getReference(next_three,3),true) ); //add each letter
        next_three += 3;
    }

    //quint64 leftover = bpos+valen - next_three;
//    cerr << "Leftover for insertion:" << leftover << endl;

    return "p."+ref_p+QString::number(position)+'_'\
            +QString::number(position+num_proteins_ins)+"ins"\
            +intermediate_proteins;

}


QString Appender::antonorakis(
       ExonData *current_exon, bool &direction, quint64 &bpos,
                              QString &VREF, QString &VALT,
                              QString &ref_codon, QString &alt_codon,
                              bool &isSNV, bool &isDel,
                              quint64 &next_codon_bpos, int regulatory)
{
    //Exon data:
    // [#Exon] -->
       // [pos1,
       //  pos2,
       //  cumulative coding length from 5->3,
       //  cumulative from 3->5]

//    cerr << direction << " " << bpos << " " << VREF.toUtf8().data() << " " << VALT.toUtf8().data() << " -- "
//         << ref_codon.toUtf8().data() << " " << alt_codon.toUtf8().data() << " - " << isSNV << " " << isDel << "   " << next_codon_bpos << endl;


//    cerr << "Exon start:" << current_exon->start << endl;
//    cerr << "Exon end:" << current_exon->end << endl;
//    cerr << "Exon 5'3:" << current_exon->five_to_three << endl;
//    cerr << "Exon 3'5:" << current_exon->three_to_five << endl;

    if (regulatory != None){
        int exon_pos = direction?current_exon->five_to_three:current_exon->three_to_five;

        int offset = -1;
        if (regulatory == Upstream){
            exon_pos -= current_exon->end - current_exon->start; // jump to start
            offset = direction?(bpos - current_exon->start):(current_exon->end - bpos);
        }
        else if (regulatory == Downstream){
            // We are at the end of exon
            offset = direction?(bpos + 1 - current_exon->end):(current_exon->start - bpos + 1);
        }

        if (direction){
            offset -= 1;
        }

        return c_nomenclature_regulatory(exon_pos, offset, VREF, VALT, isSNV, isDel)+"|| * ||* *";
    }


    int coding_pos = -1; //coding position

    if(!direction){ // reverse
        t_exon rollover = current_exon->three_to_five - (current_exon->end - current_exon->start); // #coding up to end of current exon (measuring from 3->5)
        coding_pos = (current_exon->end-bpos) + rollover + 1;  // coding from end to current
    }
    else {
        t_exon rollover = current_exon->five_to_three - (current_exon->end - current_exon->start);   // #coding up to beginning of current exon
        coding_pos = (bpos - current_exon->start) + rollover;   //coding from start to current
    }


    int protein_position = ((coding_pos-1)/3) + 1; // same for both directions (start point changes that's all, calculation same)

    //Check proteins
    QStringList proteinRA = ph->codons2Proteins(ref_codon, alt_codon).split(',');
    QString &ref_protein = proteinRA[0], &alt_protein = proteinRA[1];

    return c_nomenclature(coding_pos, VREF, VALT, isSNV, isDel)+"||"\
            +p_nomenclature(bpos, protein_position, ref_protein,alt_protein, VREF.length(), VALT.length(), isSNV, isDel, next_codon_bpos)+"||"\
            +ref_protein+' '+alt_protein;
}



#define printNewHeads {\
    if (!found_header_clist) fileout(HEADER_CLIST_FULL);\
    if (!found_header_plist) fileout(HEADER_PLIST_FULL);\
    if (!found_header_mlist) fileout(HEADER_MLIST_FULL);\
    if (!found_header_dlist) fileout(HEADER_DLIST_FULL);\
    if (!found_header_chlist) fileout(HEADER_COCHLIST_FULL);\
    if (!found_header_phlist) fileout(HEADER_PRCHLIST_FULL);\
}


void Appender::handleHeaders(QString &line, QTextStream &in, uint &countline)
{
    line="##";
    qint64 bufferpos = 0;

    bool found_header_clist = false, found_header_plist = false, found_header_mlist = false, found_header_dlist = false;
    bool found_header_chlist = false, found_header_phlist = false;
    short found_format = -1;

    bool quit_now = false;

    //Print headers
    do{
        line = in.readLine();
        //emergency brake if #CHROM not found
        if (!line.startsWith("#")) break;

        //Look for ##
        if (line.startsWith(HFORMAT)){
            found_format = 0; //Found
            if (line.startsWith(HEADER_CLIST)) found_header_clist = true;
            if (line.startsWith(HEADER_PLIST)) found_header_plist = true;
            if (line.startsWith(HEADER_MLIST)) found_header_mlist = true;
            if (line.startsWith(HEADER_DLIST)) found_header_dlist = true;
            if (line.startsWith(HEADER_COCHLIST)) found_header_chlist = true;
            if (line.startsWith(HEADER_PRCHLIST)) found_header_phlist = true;
        }
        else if (line.startsWith("#CHROM")){
            QStringList tokes = line.split('\t');
            bool form = false;
            for (int p=0; p < tokes.length(); p++){
                QString upper_header = tokes[p].toUpper();

                if (upper_header.contains("FORMAT")) {
                    FORMAT_INDEX = p;
                    INDIV_START_INDEX = p+1;
                    form = true;
                }
                else if (upper_header.contains("INFO")) {
                    INFO_INDEX = p;
                }
            }
            if (!form) {
                cerr << "\n\nCould not find FORMAT column! " << endl;
                exit(-1);
            }
            else {
                quit_now = true;
            }
        }

        else {
            if(found_format==0){ // Found format but now it's ended
                printNewHeads;  //Never found header, print new one
                found_format=-2; // so that we know that it at least exists
            }
        }
        fileout(line.toUtf8().data());
        countline++;
        bufferpos = in.pos();
    }
    while(!quit_now && !in.atEnd());

    // All headers already present? Already processed, skip.
    if (found_header_chlist && found_header_clist && found_header_dlist
            && found_header_mlist && found_header_phlist && found_header_plist){
        cerr << "Already processed, skipping." << endl;
        exit(0);
    }

    //Never found ##FORMAT in header
    if(found_format==-1){
        cerr << "No Format line in header!\nPrinting new one anyway..." << endl;
        printNewHeads;  //Never found header, print new one
    }

    //Go to end of headers
    in.seek(bufferpos);
}



bool Appender::isaSNV(QString &info_text){
    return !info_text.contains("IndelType=");
}



