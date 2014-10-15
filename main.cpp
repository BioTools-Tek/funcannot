#define VERS "v1.1"

#define G_ARG "--geneid="
#define T_ARG "--typid="

#define garg_id "AL"
#define targ_id "TYP"

#include "appender.h"

void usage(){
    cerr << "funcannot" << "  " << VERS << endl;
    cerr << "Annotates each line of a VCF file to show codon, protein, and mutation for each gene given in the genelist" << endl;
    cerr << "Please note that TYP and GENE annotations must have been performed prior to running this program" << endl;
    cerr << "\n    usage: ./funcannot <INPUTS> <FLAGS> <OUTPUTS>" << endl;
    cerr << endl;
    cerr << "(the following arguments are all MANDATORY)" << endl;
    cerr << "INPUTS:" << endl;
    cerr << "file1.vcf[+file2.vcf]  VCF file, or list seperated with '+' (NO SPACES)" << endl;
    cerr << "input.genemap          Genemap for positions of genes/exons" << endl;
    cerr << "input.dnamap           DNA codon map of format: Alu[TAB]AAC,AGC,GCA" << endl;
    cerr << "FASTA_folder           Folder containing FASTA .fa files for each chromosome" << endl;

    cerr << "\nFLAGS:" << endl;
    cerr << G_ARG << garg_id << "            specifies common genelist identifier in VCF file(s)" << endl;
    cerr << T_ARG << targ_id << "            specifies common type (SNP/Indel) identifier in VCF file(s)" << endl;

    cerr << "\nOUTPUTS:" << endl;
    cerr << "annotated_folder       each of the annotated VCF files will be placed here" << endl;
    cerr << "rejects_folder         each of the corresponding rejects will be placed here" << endl;
    cerr << endl;
    exit(-1);
}


int main(int argc, char *argv[])
{    
    if (argc<9) usage();

    QStringList vcf_files = QString(argv[1]).split('+');
    QString gmp_file = argv[2];
    QString dna_file = argv[3];
    QString fas_folder = argv[4];

    QString opt_G = QString(argv[5]).split('=')[0].trimmed()+'=';
    QString opt_T = QString(argv[6]).split('=')[0].trimmed()+'=';

    if (opt_G!=G_ARG){cerr << "should be '" << G_ARG << "' here" << endl; exit(-1);}
    if (opt_T!=T_ARG){cerr << "should be '" << T_ARG << "' here" << endl; exit(-1);}

    QString G_id = QString(argv[5]).split(G_ARG)[1].trimmed();
    QString T_id = QString(argv[6]).split(T_ARG)[1].trimmed();

    QString output_fold=argv[7];
    QString rejects_fold=argv[8];

    //Attempt to make folders -- clearly works only on unix
    QString command = "mkdir -p "+output_fold+"; mkdir -p "+rejects_fold;
    system(command.toUtf8().data());

    GeneMap gm(gmp_file);
    ProteinHandler ph(dna_file);

    for (int f=0; f< vcf_files.length(); f++){
        Appender(vcf_files[f], fas_folder, G_id, T_id, output_fold, rejects_fold, (f+1), gm, ph);
    }
    cerr << endl;
}
