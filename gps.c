#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include "nrutil.h"
#include <math.h>


// Defines
#define MAX_LATTICE 10
// Change the next three define parameters together
#define MAX_SEQ 26
#define NUMBER_OF_SPECIES 5
#define ISASSYMETRIC 0 // Not used at this point in time

#define SEQ_PER_SPECIES 5
#define MAX_NODES 2*MAX_SEQ-1
#define MAX_SEQ_LENGTH 600 // full barcode
//#define MAX_SEQ_LENGTH 150 // mini-barcode
#define MAX_TITLE_LENGTH 200
#define ITERATES 1000

struct {
    int    deme_x, deme_y;      // deme (x,y) coordinates in demes number
    double gps_x, gps_y;        // deme (w,z) coordinates in GPS hard values
    int	   lattice;		// current lattice block
    char   title[MAX_TITLE_LENGTH];
    char   simulated[MAX_SEQ_LENGTH];
    int    anc;                 // ancestor number
    int    desc1;               // descendant (from simulate.c)
    int    desc2;               // descendant (from simulate.c)
    double anc_time;            // time to its' ancestor (from simulate.c)
    double time;                // beginning time to this sample (from simulate.c)
    int    mutations;           // max possible no mutations applied (not =~ segregating sites)
    int    l_visited;
    int    r_visited;
    // double time; // time to that ancestor; ori semantic; DISCREPANCY
} sample[ MAX_NODES ];

int k_seq[MAX_NODES][MAX_NODES];                      // holds remaining seq(s) on a_lattice to undergo events (coal or mig)
int nn[MAX_NODES];                                    // counter of remaining seq(s), j, on a_lattice to undergo events
int k_lattice[MAX_NODES][MAX_LATTICE*MAX_LATTICE];    // MAX_NODES lattices, each of size MAX_LATTICE^2; physical location only
int n_lattice=0;                                      // counter of number of lattices
int k[MAX_NODES];                                     // holds an array of interior or exterior or both nodes

struct {
    int left;
    int right;
    double distance;
    int taxa_id;
    int checked;
} nj_tree[MAX_NODES];

typedef struct {
    int desc1, desc2;
} a_desc;

// Global variables
int    number_of_seq, length;
int    max_anc = 0;
int    old_anc = 0;
double alpha;                   // distance_between_nodes
long   *iseed, idum;            // For random number generator
a_desc desc;
char   DNA[4];

// for tree representation
char   newickTree[2000];
int    newickCount = 0;
int    newickRoot;
int    trackAnc[MAX_NODES];
int    trackAncCount = 0;
char   tree_output_file[256];
char   tree_accepted_output_file[256];

// external function declarations
extern float ran1(long *idum);
extern float expdev(long *idum);
extern float poidev(float xm, long *idum);

// function declarations
void   usage(void);
void   sample_coalescent(double alpha, int demes, double theta, double mig, double scale);
void   mutate_coalescent(double alpha, int demes, double theta, double mig);
void   set_scaled_positions(int begin_nn, int end_nn, int a_lattice, int demes, double scale, int sampling_scheme);
void   set_uncoalesced_position(int a_seq, int a_lattice, int demes, double scale);
int    choose_a_gene(int in_deme, int nseq, int demes, int nn, int a_lattice);
a_desc find_descendants(int a_seq);
void   inorder(int n);
void   inorder_left(int n);
void   inorder_right(int n);
void   inorder_up(int n);
void   print_block(int begin_nn, int end_nn, int a_lattice, int d, int nn); // print coalescent result of lattice
void   print_lattice(int a_lattice, int d);     // print number of and which seqs reside in each deme of the lattice
void   print_transform(int d); // transform dxd or (x,y) coordinates to k
void   print_sample(int a_seq, int demes); // print info for sequence, a_seq
void   print_mutations(int a_seq); // print history of mutations from a_seq to its' highest ancestor
// external descendants are seqs from 0 ... MAX_SEQ-1
int    find_first(int a_seq);   // find largest external descendant of a_seq
int    find_last(int a_taxon, int a_seq); // find least common ancestor of a_taxon that incl seqs (a_seq+SEQ_PER_SPECIES)
void   find_all(int a_seq);     // find all external descendants of a_seq; place in global var k[]
int    find_all_taxa1(int a_seq); // find all descendants that coal with taxa0 (seq 1-5) (a_seq = 1)
int    find_all_incl(int a_seq); // same as find_all but include internal nodes

// from simulate.c
void   do_block_coalescent(int begin_nn, int end_nn, int a_lattice, double begin_time, double end_time, int demes, double mig, double scale);

FILE   *out_debug;      // must span various functions
int    print2File = 0;

/**************************/


int main(int argc, char *argv[]) {

    /* file pointers */
    FILE   *out_log;
    FILE   *out_data;

    /* file names */
    char   filename[10][256];
    char   output_file[256];

    int    be_verbose;                          // output program status while running
    int    nflag;                               // treat N's as '-''s
    int    isSS;                                // flag to indicate ss ,
    //int    sampling_scheme;                     // sampling scheme = {0,15}
    int    d;                                   // number of demes: d*d
    double theta, mig, scale;
    int    i, j, l, m, n, anc_t, curr_t;        // counters and place holders
    int    count;                               // number of program arguments
    int    iter;                                // number of data sets to generate
    int    temp1;
    double temp2;
    int    startSeq;
    int    LCA;
    int    first,last,maxSeq;                   // simulate.c; diagnostics of query
    int    noOtherSeq;                          // number of seq from other taxa (re: taxa0)
    time_t tp;

    int printSeq = 1;                           // flag to print sample and query sequences
    int printSeq2File = 1;                      // flag to print debug output to file
    char   command[1000];                       // holds system command
    char   buffer[100];                         // holds add-on to system command
    
    int number_of_species = NUMBER_OF_SPECIES;  // for cal segregating sites
    int seg_sites[number_of_species];           // for cal segregating sites
    int seg_sites_q[number_of_species];         // for cal segregating sites
    int seg_sites_tot[number_of_species];       // for cal segregating sites

    double heights[number_of_species];
    double lengths[number_of_species];
    int    mutns[number_of_species];
    int    seq_mutns[MAX_NODES];
    double diff;
    double sum;
    double e_ss,e_height,e_length = 0.0;  // avg ss (check mutate), avg height and length (check coal)
    double avg_mutns,avg_ss,avg_height,avg_length = 0.0;
    

    
    // Model is d x d demes arranged in a lattice.  Have theta_i in
    // each deme and have migration at rate m between each deme
    // (island model = every deme is equivalent).
    // Fit the GPS coordinates onto the center of the demes with
    // maxGPS coordinates scaled to be s x d across (0<s<1; where 0
    // implies all samples from the central deme, 1 implies samples
    // scattered maximally across the demes.

    // Have theta and migration value
    // 5 sequences per group; each group occupies own lattice
    // Have time value
    // Let sequences coalesce (in each lattice?)
    // Uncoalesced sequences taken and randomly redistributed on a new
    // lattice
    // Continue until all sequences have coalesced

    d=5;                             // set as default but allow an option to change it
    be_verbose = 1;                  // set as default but allow an option to change it
    scale=0.1;                       // set as default but allow an option to change it
    number_of_seq = MAX_SEQ; // 51
    length = MAX_SEQ_LENGTH; // 600
    
    strcpy(output_file,"status_output"); // stdout output

    if(argc < 9) { usage(); exit(0);} // gps -a alpha -t theta -m migration -d demes
    for(count=1; count<argc; count++) {
        if(argv[count][0] != '-') { usage(); exit(0);}
        else switch(argv[count][1]) {
            case 'a'://alpha; time to speciation or coalescence
                alpha = atof(argv[++count]);
                continue;
            case 'd': //Number of demes or subpopulations
		d = atoi(argv[++count]);
                continue;
            case 'h': //Help
                usage();
                exit(0);
            case 'l': //GPS coordinates file
                strcpy(filename[1], argv[++count]);
                continue;
            case 'm': //Migration
                mig = atof(argv[++count]);
                continue;
//             case 'n':
//		nflag=1;// treat N's as '-''s
//                continue;
//            case 'o': // output file name
//                strcpy(filename[2], argv[++count]);
//                continue;
//            case 's': //seq file of taxa
//                strcpy(filename[0], argv[++count]);
//                continue;
            case 't':
                theta = atof(argv[++count]);
		theta = theta/MAX_SEQ_LENGTH; // per bp
                continue;
            case 'v':
                be_verbose = 0;  // turn verbosity off
                continue;
            default:
                usage();
                exit(0);
        }
    }

/*******************************************************************/
    /* Prepare the output file */
    if((out_log=fopen(output_file,"w"))==NULL) {
        fprintf(stderr,"\n Unable to write to output file %s\n\n", output_file);
        usage();
        exit(1);
    }
    if((out_data=fopen("data","w"))==NULL) {
        fprintf(stderr,"\n Unable to write to output file data\n\n");
        usage();
        exit(1);
    }
/******************************* Finished input **************************************/

    iseed = &idum;
    idum = (unsigned) time(&tp);
    idum = -idum;
    //idum = -12345; // debug purposes only
    
    //*********** Section 1: Start at some simple values
    
    //d=2;
    //theta=0.004;
    //mig = 5.0;
    
    // Output to file
    fprintf(out_log,"Setting parameters ...\n");
    fprintf(out_log,"\tSetting d to %d\n\n",d);
    fprintf(out_log,"\tSetting theta (4Nmu) to %10.5f\n\n",theta); // 4Nmu
    fprintf(out_log,"\tSetting migration (4Nm) to %10.5f\n\n",mig); // 4Nm
    fprintf(out_log,"\tSetting alpha (time to coalesce) to %10.5f\n\n",alpha);
    fprintf(out_log,"\tSetting scale (NOT used) to %10.5f\n\n",scale);
    
    fprintf(out_log,"Query is seq0 and taxa0,\n");
    i=1; j=SEQ_PER_SPECIES;
    fprintf(out_log,"First species  is taxa0 with sequences numbered from %2d to %2d\n",i,j);
    i=j+1; j=j+SEQ_PER_SPECIES;
    fprintf(out_log,"Second species is taxa1 with sequences numbered from %2d to %2d\n\n",i,j);

    fprintf(out_log,"Setting output files ...\n");
    fprintf(out_log,"\tSending status output to file: %s\n\n",output_file);

    //*********** Section 2: Generate a sample with these values

    // check avg ss, avg height and avg length
    for(i=0; i<number_of_species; i++){
        seg_sites_tot[i] = 0;
        heights[i]=0.0;
        lengths[i]=0.0;
        mutns[i] = 0;
    }
    for(i=0; i<MAX_NODES; i++){
        seq_mutns[i] = 0;
    }

    for(iter=0; iter<1; iter++){
        // initialize
        for(i=0; i<MAX_NODES; i++) {
            sample[i].desc1 = sample[i].desc2 = sample[i].anc = 0;
    	    sample[i].anc_time = 0.0;
    	    sample[i].time = 0.0;
            sample[i].lattice = -1;
            sample[i].mutations = -1;
            sample[i].l_visited = -1;
            sample[i].r_visited = -1;
        }
        for(i=0; i<MAX_NODES; i++){
            // no seqs on any of the lattices
    	    for(j=0; j<MAX_LATTICE*MAX_LATTICE; j++){
    	        k_lattice[i][j] = 0; //k_lattice[MAX_NODES][MAX_LATTICE*MAX_LATTICE]
    	    }
            for(j=0; j<MAX_NODES; j++){
                k_seq[i][j] = 0; // k_seq[MAX_NODES][MAX_NODES]
            }
            nn[i] = 0; // nn[MAX_NODES]
        }
        n_lattice = 0;
        for(i=0; i<MAX_NODES; i++){ k[i] = -1; } // clear
    
        if(d==1){
            sample_coalescent(alpha,d,theta,0.0,scale);
            mutate_coalescent(alpha,d,theta,0.0);
        } else {
            if(print2File && (out_debug=fopen("debug2file","w+"))==NULL) {
                fprintf(stderr,"\n Unable to write to output file debug2file\n");
                exit(1);
            }

            sample_coalescent(alpha,d,theta,mig,scale);
            
            //number_of_seq = MAX_SEQ;
            //length = MAX_SEQ_LENGTH;
            //if(0){ printf("max_anc: %d\tmax: %d\n",max_anc,number_of_seq*2-2); }
            
            mutate_coalescent(alpha,d,theta,mig);

            /********************************************************************/

//            for(j=1;j<MAX_SEQ;j++){
//                print_mutations(j);
//            }

            /********************************************************************/
            // calculate segregating sites for each species:
            // without query (sample[0]): seg_sites[]
            // with query: seg_sites_q[]
            for(i=0,j=1; i<number_of_species; i++){
                seg_sites[i]=0;
                seg_sites_q[i]=0;
                l=j; // start of current taxon, 1
                j=j+SEQ_PER_SPECIES; // start of next taxon, 6
                for(m=0; m < length; m++){
                    if(SEQ_PER_SPECIES > 1){ // if more than 1 seq in species
                        isSS = 0;
                        for(n=1; n<SEQ_PER_SPECIES;n++){ // for each seq in species
                            if(sample[l].simulated[m] != sample[l+n].simulated[m]
                                    && sample[l].simulated[m] != '-'
                                    && sample[l+n].simulated[m] != '-'){
                                seg_sites[i]++;
                                seg_sites_q[i]++;
                                isSS = 1;
                                //printf("@%d: %c %c, l = %d and l+j = %d: For species group %d, calculated ss=%d, ss+query=%d\n", m, sample[l].simulated[m], sample[l+n].simulated[m], l, l+n, i, seg_sites[i], seg_sites_q[i]);
                                break;
                            }
                        }// end for
                        if(isSS == 0){ // no w/i ss; check vs query
                            for(n=1; n<SEQ_PER_SPECIES;n++){ // for each seq in species
                                if(sample[l].simulated[m] != sample[0].simulated[m] // compare vs query
                                        && sample[l].simulated[m] != '-'
                                        && sample[0].simulated[m] != '-'){
                                    seg_sites_q[i]++;
                                    //printf("@%d w Query: %c %c, l = %d: For species group %d, calculated ss=%d, ss+query=%d\n", m, sample[l].simulated[m], sample[0].simulated[m], l, i, seg_sites[i], seg_sites_q[i]);
                                    break;
                                }
                            }// end for
                        }//end if
                    }else{ // 1-seq species
                        if (sample[l].simulated[m] != sample[0].simulated[m] // compare vs query
                                        && sample[l].simulated[m] != '-'
                                        && sample[0].simulated[m] != '-'){
                            seg_sites_q[i]++;
                            //printf("@%d w Query: %c %c, l = %d: For species group %d, calculated ss=%d, ss+query=%d\n", m, sample[l].simulated[m], sample[0].simulated[m], l, i, seg_sites[i], seg_sites_q[i]);
                        }
                    }
                }
            }

            /* sanity check */
            for(i=0; i<number_of_species; i++){
                if(seg_sites_q[i] > length){
                    fprintf(stderr,"Have more segregating sites than sites.\n");
                    fprintf(stderr,"Have calculated ss=%d for species %d.\n\n",seg_sites_q[i],i);
                    exit(0);
                }
                if(seg_sites_q[i] < seg_sites[i]){
                    fprintf(stderr,"Have more segregating sites without query than \n");
                    fprintf(stderr,"with it.  Not possible; there is an error somewhere.\n");
                    fprintf(stderr,"With %d; without %d for species %d.\n\n",seg_sites_q[i],seg_sites[i],i);
                    exit(0);
                }
            }
            // tally total segregating sites per species
            for(i=0; i<number_of_species; i++){ 
                seg_sites_tot[i] += seg_sites[i];
                //printf ("[%d]: %5d%5d%5d\n",i,seg_sites[i],seg_sites_q[i],seg_sites_tot[i]);
            }

            for(i=0, j=1, diff=0.0; i<number_of_species; i++){
                
                // calculate observed height; expected is E(H_n) = 2(1-1/n)
                // check tree height of each taxa
                // find LCA of taxon
                // LCAs' time = height
                LCA = find_last(i,j);
                heights[i] += sample[LCA].time;
                //printf("LCA of taxa %d (start seq: %d): %d\n",i,j,LCA);
                j+=SEQ_PER_SPECIES;

                startSeq = i*SEQ_PER_SPECIES+1; // first seq in taxon, i
                
                // calculate observed length; expected is E(L_n) = 2 sum_j=1^(n-1) 1/j
                // sum the length of all branches (only seqs in taxon)
                // for each sequence of taxon, i, get its ancestor
                // if the branch between anc and curr seq has not been visited
                // include the branch length (time difference)
                // and record the visit
                for(l=0; l<SEQ_PER_SPECIES; l++){
                    curr_t = startSeq+l;
                    // sum lengths from seq to LCA
                    do{
                        anc_t = sample[curr_t].anc; // get the ancestor
                        if(sample[anc_t].desc1 == curr_t && sample[anc_t].l_visited != 1){
                            sample[anc_t].l_visited = 1;
                            diff += sample[anc_t].time - sample[curr_t].time; // anc - current time
                        }else if(sample[anc_t].desc2 == curr_t && sample[anc_t].r_visited != 1){
                            sample[anc_t].r_visited = 1;
                            diff += sample[anc_t].time - sample[curr_t].time;
                        }else{ /* already added; do nothing */ }
                        curr_t = anc_t; // set curr equal to anc
                    }while(curr_t != LCA);
                }
                lengths[i] += diff;
                diff = 0.0;

                // calculate avg number of possible mutations
                // find LCA of taxon
                // for each sequence of taxon, i, add its mutation
                // for its ancestor (if not included), add its
                // mutation and record the visit
                // do this until LCA reached
                // NOT SURE if this works
                for(l=0; l<MAX_NODES; l++){ sample[l].l_visited = -1; sample[l].r_visited = -1; }
                for(l=0; l<SEQ_PER_SPECIES; l++){
                    curr_t = startSeq+l;
                    mutns[i] += sample[curr_t].mutations; // record mutns of curr seq
                    sample[curr_t].l_visited = 1;
                    sample[curr_t].r_visited = 1;
                    do{
                        anc_t = sample[curr_t].anc;
                        // if not counted, record mutn of anc
                        if(sample[anc_t].l_visited != -1 && sample[anc_t].r_visited != -1){
                            mutns[i] += sample[anc_t].mutations;
                            sample[anc_t].l_visited = 1;
                            sample[anc_t].r_visited = 1;
                        }
                        curr_t = anc_t;
                    }while(curr_t != LCA);
                }
            }
            // sum mutations per seq
            for(i=0; i<MAX_SEQ; i++){ seq_mutns[i] += sample[i].mutations; }
            /********************************************************************/

            // print sample sequences
            if(printSeq){
                for(i=1; i<MAX_SEQ; i++){
                    if(printSeq2File){
                        fprintf(out_data,"> %d | Genus taxa%d | 12345 | nothing \n",i, (int) ((i-0.1)/5));
                    }else{
                        fprintf(stdout,"> %d | Genus taxa%d | 12345 | nothing \n",i, (int) ((i-0.1)/5));
                    }
                    for(j=0; j<length; j++) {
                        if(printSeq2File){ fprintf(out_data,"%c",sample[i].simulated[j]);
                        }else{ fprintf(stdout,"%c",sample[i].simulated[j]); }
                    }
                    if(printSeq2File){ fprintf(out_data,"\n"); }
                    else{ fprintf(stdout,"\n"); }
                }
            }
            //print_sample(0,d);
            first = find_first(0); // find largest external desc of (anc of 0)
            last = find_last(0,0); // k[] holds all descendants of LCA of 0-5 (and any other seqs)
            for(j=1; j<MAX_SEQ; j++){ // find largest value of k[]; found at pos 0
                if (k[0]<k[j]){
                    temp1 = k[j];
                    k[j] = k[0];
                    k[0] = temp1;
                }
            }
            maxSeq = k[0];
            // - first is the largest seq that (first) coalesces with
            // seq 0 (the Query)
            // - last is the largest common ancestor (LCA) of taxa0 (6
            // seqs: 0-5)
            // - k[0] holds the largest seq that is part of coalescent
            // rooted by last (LCA); this should indicate the the
            // maximum sequence, and species, included in a coalescent
            // with species 1 (that is, identify if species 1 is
            // monophyletic)
            if(printSeq){
                if(printSeq2File){
                    fprintf(out_data,"> Query | from taxa0, 1st coalesce %2d, last coalesce with id %2d which includes max id %2d",first,last,k[0]);
                }else{
                    fprintf(stdout,"> Query | from taxa0, 1st coalesce %2d, last coalesce with id %2d which includes max id %2d",first,last,k[0]);
                }
            }

            // - first (new) is the LCA of taxa0 (5 seqs: 1-5)
            for(j=0; j<MAX_NODES; j++) { k[j]=-1; } // reset k
            first = find_all_taxa1(1);
            find_all(first); // k[] holds all descendants of LCA of 1-5 (and any other seqs)

            for(noOtherSeq=0,j=0; j<MAX_SEQ; j++){ //find all from other taxa
                // exclude the query (0)
                if(k[j] >= 6 && k[j] != 0){
                    noOtherSeq++;
//                    printf("k[%d]: %d\n", j, k[j]);
                }
            }
            if(printSeq){
                if(printSeq2File){ fprintf(out_data," - and %2d from other taxa excluding seq0\n", noOtherSeq); }
                else{ fprintf(stdout," - and %2d from other taxa excluding seq0\n", noOtherSeq); }
                for(l=0; l<length; l++){
                    if(printSeq2File){ fprintf(out_data,"%c",sample[0].simulated[l]); }
                    else{ fprintf(stdout,"%c",sample[0].simulated[l]); }
                }
                if(printSeq2File){ fprintf(out_data,"\n"); }
                else{ fprintf(stdout,"\n"); }
            }
            //if(print2File){ fclose(out_debug); } // use if exit before end of iter
            
            // print debug2file logs
            // taxa2
//            sprintf(buffer,"%d",iter);
//            if(maxSeq == 11 || maxSeq == 12 || maxSeq == 13 || maxSeq == 14 || maxSeq == 15){
//                strcat(command,"cp debug2file debug2file.taxa2.");
//                strcat(command,buffer);
//                system(command);
//            }// taxa3
//            else if(maxSeq == 16 || maxSeq == 17 || maxSeq == 18 || maxSeq == 19 || maxSeq == 20){
//                strcat(command,"cp debug2file debug2file.taxa3.");
//                strcat(command,buffer);
//                system(command);
//            }else{
//                strcat(command,"cp debug2file debug2file.iter.");
//                strcat(command,buffer);
//                system(command);
//            }
//            memset(command, '\0', sizeof(command));
//            memset(buffer, '\0', sizeof(buffer));

        }
    }// end iter<10000 or samples to generate

    // cal avg segregating sites, height and length
    // expect E(SS) = theta * sum_j=1^(n-1) 1/j
    // expect E(H_n) = 2(1-1/n)
    // expect E(L_n) = 2 sum_j=1^(n-1) 1/j
    for(i=1,sum=0.0; i<SEQ_PER_SPECIES; i++){
        sum += 1.0/i;
    }
    e_ss = (theta*length) * sum;
    e_height = 2 * (1-(1.0/SEQ_PER_SPECIES));
    e_length = 2 * sum;
    fprintf(out_log,"Data diagnostics ...\n");
    fprintf(out_log,"\ttheta * len: %5.5f *%5d = %10.5f\n",theta,length,theta*length); // 4Nmu
    fprintf(out_log,"\tsum_j=1^{n-1}; n=%d: %10.5f\n", SEQ_PER_SPECIES, sum);
    fprintf(out_log,"\tsamples: %d\n",iter);
    //printf("[E]%10s%10s%10s%10s\n","[M]","[SS]","[H]","[L]");
    //printf("[E]%20.5f%10.5f%10.5f\n",e_ss,e_height,e_length);
    fprintf(out_log,"\t[E]%10s%10s%10s\n","[SS]","[H]","[L]");
    fprintf(out_log,"\t[E]%10.5f%10.5f%10.5f\n",e_ss,e_height,e_length);
    fprintf(out_log,"\t============================================\n");
    for(i=0; i<number_of_species; i++){
        //avg_mutns = mutns[i]/(iter+0.0);
        avg_ss = seg_sites_tot[i]/(iter+0.0);
        avg_height = heights[i]/(iter+0.0);
        avg_length = lengths[i]/(iter+0.0);
        //printf("[%d]%10.5f%10.5f%10.5f%10.5f\n",i,avg_mutns,avg_ss,avg_height,avg_length);
        fprintf(out_log,"\t[%d]%10.5f%10.5f%10.5f\n",i,avg_ss,avg_height,avg_length);
    }
    /*
    printf("\nAvg mutations per sequence\n");
    printf("============================================\n");
    for(i=0; i<MAX_SEQ; i++){
        temp2 = seq_mutns[i]/(iter+1.0);
        printf("[%d]:%10.5f\n",i,temp2);
    }
    */
    fclose(out_log);
    fclose(out_data);
    if(print2File){ fclose(out_debug); }
    return 0;
}//end main


/********************************************************************/
void usage(void) {
    fprintf(stderr,"Syntax: gps -d deme -t theta -m migration -a alpha\n");
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"\t -d #; number of demes [default = 5] \n");
    fprintf(stderr,"\t -t theta (per locus) mutation rate;\n");
    fprintf(stderr,"\t -m migration (4Nm);\n");
    fprintf(stderr,"\t -a distance_between_nodes; time to speciation/coalescence;\n");
//    fprintf(stderr,"\t -s taxa_input.filename;\n");
//    fprintf(stderr,"\t\tformat line 1: # of taxa and # length of seq \n");
//    fprintf(stderr,"\t\tformat line >=2: title[10chars]/space/sequence all on one line \n");
//    fprintf(stderr,"\t -n; treat N's in the sequence as '-'s\n");
    fprintf(stderr,"\t -l gps_data.filename; (irrelevant at the moment; pairs of unused numbers) \n");
//    fprintf(stderr,"\t -o filename; output file [default 'output']\n");
//    fprintf(stderr,"\t -v; turn off verbosity [default on]\n");
    fprintf(stderr,"\t -h; print this information\n");
    fprintf(stderr,"\n");
}
/********************************************************************/
void  sample_coalescent(double alpha, int d, double theta, double mig, double scale) {
    int debug = 1; // DEBUG
    int sampling_scheme = 0;
    max_anc = number_of_seq;
    /*
    sampling_scheme:
        0: random (default)
        1: all
        2: 1other
        3: 2other
        4: wQ
        5: Qcloser_1other
        6: Qcloser_2other
        7: all_large
        8: 1other_large
        9: 2other_large
        10: 1other_large_rand
        11: 2other_large_rand
        12: 3other_large_rand
        13: 4other_large_rand
        14: 5other_large_rand
        15: all_large_rand
        16: gps_all_large_4x4
        17: gps_1other_large_4x4
    */
    if(ISASSYMETRIC==1){
	    sampling_scheme = 0; // see above
	    set_scaled_positions(0,5,n_lattice,d,scale,sampling_scheme);
	    if(debug){ print_lattice(n_lattice,d); }
	    do_block_coalescent(0,5,n_lattice,0.0,alpha,d,mig,scale); //0
	    if(debug){ print_block(0,5,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    // for other 9 species, place sequences randomly on each lattice
	    sampling_scheme = 0;
	    set_scaled_positions(6,10,n_lattice,d,scale,sampling_scheme);
	    if(debug){ print_lattice(n_lattice,d); }
	    do_block_coalescent(6,10,n_lattice,0.0,alpha,d,mig,scale); //1
	    if(debug){ print_block(6,10,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,10,n_lattice,alpha,2.0*alpha,d,mig,scale); //2
	    if(debug){ print_lattice(n_lattice,d); }
	    if(debug){ print_block(0,10,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(11,15,n_lattice,d,scale,sampling_scheme);
	    if(debug){
	        print_lattice(n_lattice,d);
	    }
	    do_block_coalescent(11,15,n_lattice,0.0,2.0*alpha,d,mig,scale); //3
	    if(debug){ print_block(11,15,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,15,n_lattice,2.0*alpha,3.0*alpha,d,mig,scale); //4
	    if(debug){ print_block(0,15,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(16,20,n_lattice,d,scale,sampling_scheme); //5
	    do_block_coalescent(16,20,n_lattice,0.0,3.0*alpha,d,mig,scale);
	    if(debug){ print_block(16,20,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,20,n_lattice,3.0*alpha,4.0*alpha,d,mig,scale); //6
	    if(debug){ print_block(0,20,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(21,25,n_lattice,d,scale,sampling_scheme);
	    do_block_coalescent(21,25,n_lattice,0.0,4.0*alpha,d,mig,scale); //7
	    if(debug){ print_block(21,25,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,25,n_lattice,4.0*alpha,5.0*alpha,d,mig,scale); //8
	    if(debug){ print_block(0,25,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(26,30,n_lattice,d,scale,sampling_scheme);
	    do_block_coalescent(26,30,n_lattice,0.0,5.0*alpha,d,mig,scale); //9
	    if(debug){ print_block(26,30,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,30,n_lattice,5.0*alpha,6.0*alpha,d,mig,scale); //10
	    if(debug){ print_block(0,30,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(31,35,n_lattice,d,scale,sampling_scheme);
	    do_block_coalescent(31,35,n_lattice,0.0,6.0*alpha,d,mig,scale); //11
	    if(debug){ print_block(31,35,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,35,n_lattice,6.0*alpha,7.0*alpha,d,mig,scale); //12
	    if(debug){ print_block(0,30,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(36,40,n_lattice,d,scale,sampling_scheme);
	    do_block_coalescent(36,40,n_lattice,0.0,7.0*alpha,d,mig,scale); //13
	    if(debug){ print_block(36,40,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,40,n_lattice,7.0*alpha,8.0*alpha,d,mig,scale); //14
	    if(debug){ print_block(0,40,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(41,45,n_lattice,d,scale,sampling_scheme);
	    do_block_coalescent(41,45,n_lattice,0.0,8.0*alpha,d,mig,scale); //15
	    if(debug){ print_block(41,45,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	    
	    do_block_coalescent(0,45,n_lattice,8.0*alpha,9.0*alpha,d,mig,scale); //16
	    if(debug){ print_block(0,45,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    set_scaled_positions(46,50,n_lattice,d,scale,sampling_scheme);
	    do_block_coalescent(46,50,n_lattice,0.0,9.0*alpha,d,mig,scale); //17
	    if(debug){ print_block(46,50,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    do_block_coalescent(0,50,n_lattice,9.0*alpha,100000.0*alpha,d,mig,scale); //18
	    if(debug){ print_block(0,50,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
    }else{
        sampling_scheme = 0; // see above
	    set_scaled_positions(0,5,n_lattice,d,scale,sampling_scheme);
	    if(debug){ print_lattice(n_lattice,d); }
	    do_block_coalescent(0,5,n_lattice,0.0,alpha,d,mig,scale); //0
	    if(debug){ print_block(0,5,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
	    // for other 4 species, place sequences randomly on each lattice
	    sampling_scheme = 0;
	    set_scaled_positions(6,10,n_lattice,d,scale,sampling_scheme);
	    if(debug){ print_lattice(n_lattice,d); }
	    do_block_coalescent(6,10,n_lattice,0.0,alpha,d,mig,scale); //1
	    if(debug){ print_block(6,10,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
//	    do_block_coalescent(0,10,n_lattice,alpha,2.0*alpha,d,mig,scale); //2
//	    if(debug){ print_lattice(n_lattice,d); }
//	    if(debug){ print_block(0,10,n_lattice,d,nn[n_lattice]); }
//	    n_lattice++;
	
        set_scaled_positions(11,15,n_lattice,d,scale,sampling_scheme);
	    if(debug){ print_lattice(n_lattice,d); }
	    do_block_coalescent(11,15,n_lattice,0.0,alpha,d,mig,scale); //3
	    if(debug){ print_block(11,15,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
//	    do_block_coalescent(0,15,n_lattice,2.0*alpha,3.0*alpha,d,mig,scale); //4
//	    if(debug){ print_block(0,15,n_lattice,d,nn[n_lattice]); }
//	    n_lattice++;
	
	    set_scaled_positions(16,20,n_lattice,d,scale,sampling_scheme); //5
        if(debug){ print_lattice(n_lattice,d); }
	    do_block_coalescent(16,20,n_lattice,0.0,alpha,d,mig,scale);
	    if(debug){ print_block(16,20,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
	
    	set_scaled_positions(21,25,n_lattice,d,scale,sampling_scheme);
        if(debug){ print_lattice(n_lattice,d); }
	    do_block_coalescent(21,25,n_lattice,0.0,alpha,d,mig,scale); //6
	    if(debug){ print_block(21,25,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;

	    do_block_coalescent(0,25,n_lattice,alpha,100000.0*alpha,d,mig,scale); //7
        if(debug){ print_lattice(n_lattice,d); }
	    if(debug){ print_block(0,25,n_lattice,d,nn[n_lattice]); }
	    n_lattice++;
    }
}
/********************************************************************/
void   do_block_coalescent(int begin_nn, int end_nn, int a_lattice, double begin_time, double end_time, int d, double mig, double scale) {
    int i, j, x, y;
    int c1, c2, id;
    double Ic, Im, sum_k_1, sum_k;
    double a_random, sum_migr, sum_coal;
    float event_time;
    int debug = 0;
    
    old_anc = max_anc;
    
    nn[a_lattice]=0;
    for(i=begin_nn; i<=end_nn; i++) {
	if(sample[i].anc_time == begin_time){ // (leaf) seqs that did not previously coalesce
            set_uncoalesced_position(i,a_lattice,d,scale);
            k_seq[a_lattice][nn[a_lattice]]=i; // k_seq list holds seqs currently on lattice
            nn[a_lattice]=nn[a_lattice]+1;
        }
    }
    for(i=MAX_SEQ; i<max_anc; i++) {
	if(sample[i].anc_time == begin_time){ // internal nodes that did not previously coalesce
            set_uncoalesced_position(i,a_lattice,d,scale);
            k_seq[a_lattice][nn[a_lattice]]=i;
            nn[a_lattice]=nn[a_lattice]+1;
        }
    }
    // should now have nn[a_lattice] taxa in array k_seq (holds seqs
    // currently on lattice
    
    event_time = begin_time;
    while(nn[a_lattice] > 1) {
        // Coalescent rate is of seqs in THIS lattice
	for(sum_k_1=0.0, i=0; i<d*d; i++) { sum_k_1 += 0.5*k_lattice[a_lattice][i]*(k_lattice[a_lattice][i]-1.0); }
	Ic = d*d * sum_k_1;
	// Migration rate is
	for(sum_k=0.0, i=0; i<d*d; i++) { sum_k += k_lattice[a_lattice][i]; }
	Im = 0.5*mig*d*d*sum_k;
        // time scaled in units of 2Nd*d generations
        // Ic + Im is the avg rate until (a/the first) event (units: 2Nd*d)
	event_time += expdev(iseed)/(Ic+Im);
        if(debug){ printf("event time: %.5f\tbegin time: %.5f\tend time: %.5f\n", event_time, begin_time, end_time); }

        if (event_time > end_time) { // do not coalesce or migrate if outside time bounds
            if(debug){ printf("no events\n"); }
            for(i=begin_nn; i<=end_nn; i++) { if(sample[i].anc_time == begin_time){ sample[i].anc_time = end_time; }}
            // any new ancestors that have not coalesced should be time extended
            for(i=MAX_SEQ; i<old_anc; i++) { if(sample[i].anc_time == begin_time) sample[i].anc_time = end_time; }
            for(i=old_anc; i<max_anc; i++) { if(sample[i].anc_time <= begin_time) sample[i].anc_time = end_time; }
            return;
        }
        
	if(ran1(iseed) < Ic/(Ic+Im)) { // next event is coalescence
            if(debug){ printf("coal event (%.2f) vs mig (%.2f)\n", Ic/(Ic+Im), Im/(Ic+Im)); }
	    // choose deme for event with prob kj(kj-1)/sumki(ki-1)
            // a deme is chosen weighted according to the coalescence
            // rates
	    a_random = ran1(iseed);
	    for(sum_coal=0.0, j=0; j<d*d; j++) {
                //calculate probability of coalescence for each deme
                //in lattice
		sum_coal += 0.5*(k_lattice[a_lattice][j]*(k_lattice[a_lattice][j]-1.0))/sum_k_1;
                if (debug){ printf("Pr(coal @ deme[%d]): %5.5f with %d seqs\n", j, sum_coal, k_lattice[a_lattice][j]); }
		if(a_random < sum_coal) { break; }
	    }
	    k_lattice[a_lattice][j]--;
	    nn[a_lattice]--;
            if(debug){
                for(i=0; i<d*d; i++) { printf("%d ", k_lattice[a_lattice][i]); }
	        printf("\n");
                for(i=0; i<(nn[a_lattice]+1); i++){
                    id = k_seq[a_lattice][i];
                    printf("k_seq[%d][%d]: %d\tdeme: %d\n", a_lattice, i, id, d*sample[id].deme_x + sample[id].deme_y);
                }
            }
            // set location of max_anc
	    sample[max_anc].deme_x = j/d;
	    sample[max_anc].deme_y = j%d;
            sample[max_anc].lattice = a_lattice;
	    // choose two within that deme c1, c2
	    c1 = choose_a_gene(j,k_lattice[a_lattice][j]+1,d, nn[a_lattice]+1, a_lattice);
	    do { c2 = choose_a_gene(j,k_lattice[a_lattice][j]+1,d, nn[a_lattice]+1, a_lattice); } while (c2==c1);
            if(debug){
                printf("anc id: %d\t", max_anc);
                printf("gene c1: %d\t", c1);
                printf("gene c2: %d\n\n", c2);
            }
	    sample[c1].anc = max_anc;
	    sample[c2].anc = max_anc;
	    sample[c1].anc_time = event_time;
            sample[c2].anc_time = event_time;
            // take genes c1, c2 off the grid
	    sample[c1].deme_x = sample[c1].deme_y = -1;
	    sample[c2].deme_x = sample[c2].deme_y = -1;
            for(i=j=0; j<nn[a_lattice]+1; j++) { // create a new list without c1, c2
	        if(k_seq[a_lattice][j]!=c1 && k_seq[a_lattice][j]!=c2) {
		    k_seq[a_lattice][i]=k_seq[a_lattice][j];
		    i++;
	        }
	    }
            k_seq[a_lattice][i]=max_anc;
            k_seq[a_lattice][i+1] = -1;  // just to help debug
            sample[max_anc].desc1 = c1;
            sample[max_anc].desc2 = c2;
	    sample[max_anc].time = event_time;
            max_anc++;
	} else { // next event is migration, choose a random gene to migrate
            if (debug){ printf("mig event (%.2f) vs coal (%.2f)\n", Im/(Ic+Im), Ic/(Ic+Im)); }
	    a_random = ran1(iseed);
	    for(sum_migr=0.0, j=0; j<d*d; j++) {
		sum_migr += ((double) k_lattice[a_lattice][j])/((double) nn[a_lattice]); //probability of mig
		if(a_random < sum_migr) { break; }
	    }
	    // choose one gene to move from j to new location, c1
	    c1 = choose_a_gene(j,k_lattice[a_lattice][j],d, nn[a_lattice], a_lattice);
            if(debug){ printf("choose gene @ [%d]: %d\n", j, c1); }
	    k_lattice[a_lattice][j]--;
	    // choose new i+/-1, j+/-1 location
	    // position j is position (j/d,j%d)
   	    do {
   		x = j/d;  y = j%d;
   		a_random = ran1(iseed);
   		if(a_random < 0.25) { 
   		    x++;   // i+1
   		} else { 
   		    if(a_random < 0.50) { 
   			x--;  // i-1
   		    } else { 
   			if(a_random < 0.75) { 
   			    y++; // j+1
   			} else { y--; } // j-1
   		    }
   		}
   	    } while( x >= d || x < 0 || y >= d || y < 0);
            
            //  Instead just choose it randomly.  // GBG debug
            //  but ensure is not where it started // GBG debug
//            do {
//		x = d*ran1(iseed);  // GBG debug
//		y = d*ran1(iseed);  // GBG debug
//	    } while( x==j/d && y==j%d );
            
   	    j = x*d + y; 
	    k_lattice[a_lattice][j]++;
	    sample[c1].deme_x = x; // update current location of chosen gene
	    sample[c1].deme_y = y;
            if(debug){
                for(i=0; i<d*d; i++) { printf("%d ", k_lattice[a_lattice][i]); }
	        printf("\n\n");
            }
	}// end else
    }
    sample[max_anc-1].anc_time = end_time; // have coalesced them all; only one left
    if(debug){ printf("all coalesced\n"); }
    return;
}
/********************************************************************/
void   mutate_coalescent(double alpha, int d,double theta,double mig) {
    int debug = 0;
    char X, DNA[4];
    int i,k,l;
    double time_left, time_right;
    float mean_mut_rate_left, mean_mut_rate_right;
    int mut_left, mut_right;
    int mut_left_tot, mut_right_tot = 0;

    float temp1, temp2 =0.0;
    
    DNA[0] = 'A'; DNA[1] = 'C'; DNA[2] = 'G'; DNA[3] = 'T';
    for(i=0; i<length; i++) { // generate random sequence
	l = ran1(iseed)*4;
	sample[number_of_seq*2 - 2].simulated[i] = DNA[l];
    }

    // for each new node
    for(i=number_of_seq*2 - 2; i>=number_of_seq; i--) {
	desc=find_descendants(i);
        if (debug){
            if(print2File){
                fprintf(out_debug,"ANCestor (%d) has descendants (%d) and (%d)\n", i, desc.desc1, desc.desc2);
            }else{
                printf("ANCestor (%d) has descendants (%d) and (%d)\n", i, desc.desc1, desc.desc2);
            }
        }
        
	for(k=0; k<length; k++) { // descendants a,b get copies of ancestor but it is mutated according to theta and times
	    sample[desc.desc1].simulated[k]=sample[i].simulated[k];
	    sample[desc.desc2].simulated[k]=sample[i].simulated[k];
	}
        // time b/t ancestor and each of its descendants
        // scaled by 2N(d*d) generations
	time_left = sample[i].time - sample[desc.desc1].time;
	time_right = sample[i].time - sample[desc.desc2].time;

        if (isinf(sample[i].time)){
            // special case when event_time == inf
            // generate completely new simulated sequences for
            // descendents so anc and desc do not look like each other
            // and biologically represent a situation where they
            // produce descendants that are radically different
            for(k=0; k<length; k++) {
	        l = ran1(iseed)*4;
	        sample[desc.desc1].simulated[k] = DNA[l];
	        l = ran1(iseed)*4;
	        sample[desc.desc2].simulated[k] = DNA[l];
            }
        }else{

            mut_left_tot = mut_right_tot = 0;
            
            // scaled theta ((4N_mu * d*d * k)/2) = substitutions/site/generation (scaled 2N_d*d generations)
            // scaled_theta = (theta * d*d * number_of_seq)/2; // NOT USED

            // calculate the number of.mutations mutations (poisson distributed)
            // calculate average scaled mutation rate:
            //          expected coalescence time = 2N
            //          theta == the expected no of mutations separating a sample of two sequences
            //          = 2Nu*mu (or 4Nu*mu)
            //          so, expected no of mutations for one sequence: 2Nu*mu/2 (4Nu*mu/2)
            // consider (d*d*number_of_seq) disjoint lineages: (theta * d*d * k)/2; units (2N_dd generations)
            // where theta = 4Nu*mu (pg. 109)
            
            
            // calculate mean rate of mutation: units (substitutions/gene)
            // mean rate of mutation events     = (theta (sub/bp) * length (bp/gene) * time)/2
            //                                  = (theta (sub/gene) * time)/2
            // units: substitutions/gene/2N_dd generations (is that right?) * 2N_dd generations = substitutions/gene

            // can assume poidev distribution of mutations per site
            // ONLY b/c the rate of mutation is low (expect it to be
            // comparable to when you apply number of mutations over
            // the entire length of the gene)
            // for each site
            
            /* finite sites model*/
            /*
            mean_mut_rate_left = (theta * time_left)/2;
            mean_mut_rate_right = (theta * time_right)/2;
            for(k=0; k<length; k++){

                mut_left = poidev(mean_mut_rate_left,iseed);
                mut_right = poidev(mean_mut_rate_right,iseed);
                
                for(l=0; l<mut_left;l++){
                    // mutate site of left seq mut_left times
                    do { X=DNA[(int) (ran1(iseed)*4)]; } while(X==sample[desc.desc1].simulated[k]);
	            sample[desc.desc1].simulated[k]=X;

                    //X=DNA[(int) (ran1(iseed)*4)];
                    //sample[desc.desc1].simulated[k]=X;
                }
                for(l=0; l<mut_right;l++){
                    // mutate site of right seq mut_right times
                    do { X=DNA[(int) (ran1(iseed)*4)]; } while(X==sample[desc.desc2].simulated[k]);
	            sample[desc.desc2].simulated[k]=X;

                    //X=DNA[(int) (ran1(iseed)*4)];
                    //sample[desc.desc2].simulated[k]=X;
                }
                if(mut_left != 0){ mut_left_tot++; } // no of mutations for the seq
                if(mut_right != 0){ mut_right_tot++; }
            }
            sample[desc.desc1].mutations = mut_left_tot;
            sample[desc.desc2].mutations = mut_right_tot;
            */
            /* infinite alleles model */
            mean_mut_rate_left = (theta * length * time_left)/2;
            mean_mut_rate_right = (theta * length * time_right)/2;

            mut_left = poidev(mean_mut_rate_left,iseed); // change float to int automatic
            mut_right = poidev(mean_mut_rate_right,iseed);

            sample[desc.desc1].mutations = mut_left;
            sample[desc.desc2].mutations = mut_right;
            
	    for(k=0; k<mut_left; k++) { 
                l=ran1(iseed)*length;
	        do { X=DNA[(int) (ran1(iseed)*4)]; } while(X==sample[desc.desc1].simulated[l]);
	        sample[desc.desc1].simulated[l]=X;
	    }
	    for(k=0; k<mut_right; k++) { 
                l=ran1(iseed)*length;
	        do { X=DNA[(int) (ran1(iseed)*4)]; } while(X==sample[desc.desc2].simulated[l]);
	        sample[desc.desc2].simulated[l]=X;
	    }

            if (debug){
                if(print2File){
                    fprintf(out_debug,"theta: %10.5f\n", theta);
                    fprintf(out_debug,"time_left: %10.5f\ttime_right: %10.5f\n", time_left, time_right);
                    fprintf(out_debug,"mean_mut_rate_left: %10.5f\tmean_mut_rate_right: %10.5f\n", mean_mut_rate_left, mean_mut_rate_right);
                    fprintf(out_debug,"poidev_mut_left: %d\tpoidev_mut_right: %d\n", mut_left, mut_right);
	            fprintf(out_debug,"Mutate %d randomly at %d sites\n", desc.desc1, mut_left);
	            fprintf(out_debug,"Mutate %d randomly at %d sites\n", desc.desc2, mut_right);
                }else{
                    printf("theta: %10.5f\n", theta);
                    printf("time_left: %10.5f\ttime_right: %10.5f\n", time_left, time_right);
                    printf("mean_mut_rate_left: %10.5f\tmean_mut_rate_right: %10.5f\n", mean_mut_rate_left, mean_mut_rate_right);
                    printf("poidev_mut_left: %d\tpoidev_mut_right: %d\n", mut_left, mut_right);
	            //printf("Mutate %d randomly at %d sites\n", desc.desc1, mut_left);
	            //printf("Mutate %d randomly at %d sites\n", desc.desc2, mut_right);
       	            printf("Mutate %d randomly at %d sites\n", desc.desc1, mut_left_tot);
	            printf("Mutate %d randomly at %d sites\n", desc.desc2, mut_right_tot);
                }
            }

        }//end else
    }
    for(i=number_of_seq*2 - 2; i>=0; i--) { // sanity check
        if(sample[i].simulated[0] != 'A' &&
           sample[i].simulated[0] != 'C' &&
           sample[i].simulated[0] != 'G' &&
           sample[i].simulated[0] != 'T') {
	    printf("\n");
	    printf("Something wrong with sample i=%d simulated sequence \n",i);
	    printf("%s",sample[i].simulated);
	    printf("\n");
	    printf("Stopping!\n");
	    exit(0);
	}
    }
    return;
}
/********************************************************************/
void   set_scaled_positions(int begin_nn, int end_nn, int a_lattice, int d, double scale, int sampling_scheme) {
    int i, j, l, m, n;
    int r1,r2,r3,r4,r5; // ref seq
    int debug = 0;

    /* data structure:
    int    deme_x, deme_y;      // deme (x,y) coordinates in demes number
    double gps_x, gps_y;        // deme (w,z) coordinates in GPS hard values
    
    sampling_scheme:
        0: random (default)
        1: all
        2: 1other
        3: 2other
        4: wQ
        5: Qcloser_1other
        6: Qcloser_2other
        7: all_large
        8: 1other_large
        9: 2other_large
        10: 1other_large_rand
        11: 2other_large_rand
        12: 3other_large_rand
        13: 4other_large_rand
        14: 5other_large_rand
        16: gps_all_large_4x4
        17: gps_1other_large_4x4
    */
    //-----------------------------------------------------------------------
    switch(sampling_scheme) {
        case 1:
	    // Place sequences all in one deme of the lattice (0,0) except for
	    // Q (1,1)
	    // gps_all
            if(debug){ printf ("gps all\n"); }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = sample[i].deme_y = 1; //(1,1)
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            sample[i].lattice = a_lattice;
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 2:
	    // Place most sequences all in one deme of the lattice (0,0)
	    // except for 1 randomly chosen sequence of the species (0,1) and
	    // Q (1,1)
	    // gps_1other
            if(debug){ printf ("gps 1other\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    // the chosen random seq
	    if(begin_nn == 0){ j = begin_nn + 1 + j; } // do not start from 0 (query)
	    else{ j = begin_nn + j; }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){ sample[i].deme_x = sample[i].deme_y = 1; // (1,1)
	        }else if(i==j){ //(0,1)   
	            sample[i].deme_x = 0;
	            sample[i].deme_y = 1;
	        }else{ sample[i].deme_x = sample[i].deme_y = 0; } // (0,0)
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 3:
	    // Place most sequences all in one deme of the lattice (0,0)
	    // except for 1 randomly chosen sequence in (0,1) and
	    // 1 randomly chosen (but diff) sequence in (1,0) and
	    // Q (1,1)
	    // gps_2other
            if (debug) { printf ("gps 2other\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose 1st random seq on a scale of 0-4
	    do{ l = SEQ_PER_SPECIES*ran1(iseed); // choose 2nd random seq that has not been chosen
	    }while(j==l); // ensure j != l
	    
	    // do not start from 0 (query)
	    if(begin_nn == 0){
	        j = begin_nn + 1 + j;
	        l = begin_nn + 1 + l;
	    }else{
	        j = begin_nn + j;
	        l = begin_nn + l;
	    }
	    
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){ sample[i].deme_x = sample[i].deme_y = 1; // (1,1) or 3
	        }else if(i==j){ //(0,1) or 1
	            sample[i].deme_x = 0;
	            sample[i].deme_y = 1;
	        }else if(i==l){ // (1,0) or 2
	            sample[i].deme_x = 1;
	            sample[i].deme_y = 0;
	        }else{ sample[i].deme_x = sample[i].deme_y = 0; } // (0,0) or 0
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 4:
	    /* The following scenario is the ``ultimate'' test case showing
	     * that assignment is really good when the query is located in the
	     * same place in which there is a reference sequence. This case
	     * should do better than the above three scenarios and below two
	     * scenarios. */
	    // Place most sequences all in one deme (0,0)
	    // except for 1 randomly chosen sequence to place with the Q (1,1)
	    // gps_wQ
            if(debug){ printf("gps wQ\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    // the chosen random seq
	    if(begin_nn == 0){ j = begin_nn + 1 + j; } // do not start from 0 (query)
	    else{ j = begin_nn + j; }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){ sample[i].deme_x = sample[i].deme_y = 1; // (1,1)
	        }else if(i==j){ // with Q in (1,1)
	            sample[i].deme_x = sample[i].deme_y = 1;
	        }else{ sample[i].deme_x = sample[i].deme_y = 0; } // (0,0)
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 5:
	    /* The following two scenarios tests to see if the assignment is
	     * better if the Query is located closer to most of the sample
	     * sequences.  They can be compared to scenario 3 (gps_2other). */
	    // Place most sequences all in one deme (0,0)
	    // except for 1 randomly chosen sequence (1,1)
	    // and Q is closer to most of the sample sequences (1,0) or (0,1)
	    // gps_Qcloser_1other
            if(debug){ printf("gps Qcloser 1other\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    // the chosen random seq
	    if(begin_nn == 0){ j = begin_nn + 1 + j; } // do not start from 0 (query)
	    else{ j = begin_nn + j; }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            do{ m = d*d*ran1(iseed); } // spans 0 .. 3 if d x d = 4
	            while(m==0 || m==3);
	            if(m==1){ // 2*0 + 1 = 1 or (0,1)
	                sample[i].deme_x = 0;
	                sample[i].deme_y = 1;
	            }else if(m==2){ // 2*1 + 0 = 2 or (1,0)
	                sample[i].deme_x = 1;
	                sample[i].deme_y = 0;
	            }else{ printf("ERROR: d is not 1 or 2\n"); exit(1); }
	        }else if(i==j){ // (1,1)
	            sample[i].deme_x = sample[i].deme_y = 1;
	        }else{ sample[i].deme_x = sample[i].deme_y = 0; } // (0,0)
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 6:
	    // Place most sequences all in one deme (0,0)
	    // except for 1 randomly chosen sequence (1,1) and
	    // 1 randomly chosen (but diff) sequence (0,1) or (1,0) and
	    // Q is closer to most of the sample sequences but not in the same
	    // location as the 2nd randomly chosen sequence (1,0) or (0,1)
	    // respectively
	    // gps_Qcloser_2other
            if(debug){ printf("gps Qcloser 2other\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    do{ l = SEQ_PER_SPECIES*ran1(iseed); // choose 2nd random seq that has not been chosen
	    }while(j==l); // ensure j != l
	    // do not start from 0 (query)
	    if(begin_nn == 0){
	        j = begin_nn + 1 + j;
	        l = begin_nn + 1 + l;
	    }else{
	        j = begin_nn + j;
	        l = begin_nn + l;
	    }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            do{ m = d*d*ran1(iseed); } // spans 0 .. 3
	            while(m==0 || m==3);
	            if(m==1){ // Q in (0,1)
	                sample[i].deme_x = 0; // place Q
	                sample[i].deme_y = 1;
	            }else if(m==2){ // Q in (1,0)
	                sample[i].deme_x = 1; // place Q
	                sample[i].deme_y = 0;
	            }else{ printf("ERROR: m is not 1 or 2\n"); exit(1); }
	        }else if(i==j){ // (1,1)
	            sample[i].deme_x = sample[i].deme_y = 1;
	        }else if(i==l){
	            if(m==1){ // m SHOULD be defined by now
	                sample[l].deme_x = 1; // l in (1,0)
	                sample[l].deme_y = 0;
	            }else if (m==2){
	                sample[l].deme_x = 0; // l in (0,1)
	                sample[l].deme_y = 1;    
	            }else{ // no query, randomly place second seq
	                do{ n = d*d*ran1(iseed); }
	                while(n==0 || n == 3); // want deme 1 or 2
	                if(n==1){
	                    sample[l].deme_x = 1; // l in (1,0)
	                    sample[l].deme_y = 0;
	                }else if(n==2){
	                    sample[l].deme_x = 0; // l in (0,1)
	                    sample[l].deme_y = 1;
	                }else{ printf("ERROR: should not reach here\n"); exit(1); }
	            }
	        }else{ sample[i].deme_x = sample[i].deme_y = 0; } // (0,0)
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 7:
	    /* The following tests to see what happens to the assignment if
	     * seqs are in a larger lattice.  I would guess that the
	     * assignment is no different than when d = 2 (2x2 lattice). */ 
	    // Place most sequences all in one deme (0,0) and
	    // Q in (2,0)
	    // gps_all_large
            if(debug){ printf("gps all large\n"); }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = 2;
	            sample[i].deme_y = 0;
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            sample[i].lattice = a_lattice;
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 8:
	    // Place most sequences all in one deme (0,0)
	    // except for 1 randomly chosen sequence (2,2) and
	    // Q in (2,0)
	    // gps_1other_large
            if(debug){ printf("gps 1other large\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    // the chosen random seq
	    if(begin_nn == 0){ j = begin_nn + 1 + j; } // do not start from 0 (query)
	    else{ j = begin_nn + j; }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = 2;
	            sample[i].deme_y = 0;
	        }else if(i==j){
	            sample[i].deme_x = sample[i].deme_y = 2; //(2,2)
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            sample[i].lattice = a_lattice;
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 9:
	    // Place most sequences all in one deme (0,0)
	    // except for 1 randomly chosen sequence (1,1) and
	    // 1 randomly chosen (but diff) sequence (2,2) and
	    // Q in (2,0)
	    // gps_2other_large
            if(debug){ printf("gps 2other large\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    do{ l = SEQ_PER_SPECIES*ran1(iseed); // choose 2nd random seq that has not been chosen
	    }while(j==l); // ensure j != l
	    // do not start from 0 (query)
	    if(begin_nn == 0){
	        j = begin_nn + 1 + j;
	        l = begin_nn + 1 + l;
	    }else{
	        j = begin_nn + j;
	        l = begin_nn + l;
	    }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = 2; // (2,0)
	            sample[i].deme_y = 0;
	        }else if(i==j){
	            sample[i].deme_x = sample[i].deme_y = 2; //(2,2)
	        }else if(i==l){
	            sample[i].deme_x = sample[i].deme_y = 1; //(1,1)
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            sample[i].lattice = a_lattice;
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;
        
        case 10:
	    // Place most sequences all in one deme (0,0)
	    // and place 1 random sequence in a random deme
            // and Q in (2,2)
	    // dxd = 9
	    // gps_1other_large_rand
            if(debug){ printf("gps 1other large rand\n"); }
	    r1 = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    // the chosen random seq
	    if(begin_nn == 0){ r1 = begin_nn + 1 + r1; } // do not start from 0 (query)
	    else{ r1 = begin_nn + r1; }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = sample[i].deme_y = 2; //(2,2)
	        }else if(i==r1){ // random dispersed ref seq
                // choose random deme (x,y)
                do{
                    j = d*ran1(iseed);
                    l = d*ran1(iseed);
                }while ((j == 0 && l == 0) || (j == 2 && l == 2)); // continue if (0,0) or (2,2)
                sample[i].deme_x = j;
                sample[i].deme_y = l;
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 11:
	    // Place most sequences all in one deme (0,0)
        // and 2 random sequences in random, distinct, demes
        // and Q in (2,2)
	    // d x d = 9
	    // gps_2other_large_rand
            if(debug){ printf("gps 2other large rand\n"); }
	    r1 = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    do{ r2 = SEQ_PER_SPECIES*ran1(iseed); // choose 2nd random seq that has not been chosen
	    }while(r1==r2); // ensure r1 != r2
	    // do not start from 0 (query)
	    if(begin_nn == 0){
	        r1 = begin_nn + 1 + r1;
	        r2 = begin_nn + 1 + r2;
	    }else{
	        r1 = begin_nn + r1;
	        r2 = begin_nn + r2;
	    }
        // choose random deme for each random seq
        for(i=0; i<=1; i++){
            do{ // choose random deme
                j = d*ran1(iseed);
                l = d*ran1(iseed);
                m = d*j + l; // get deme
            }while((j == 0 && l == 0) || (j == 2 && l == 2) || k_lattice[a_lattice][m] != 0);
            // choose new random deme if is a) base deme, b) deme with query and c) is occupied
            if(i==0){
                sample[r1].deme_x = j;
                sample[r1].deme_y = l;
            }else if (i==1){
                sample[r2].deme_x = j;
                sample[r2].deme_y = l;
            }
	        k_lattice[a_lattice][m]++;
        } // end random seq
        // place non-random seqs
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = sample[i].deme_y = 2; //(2,2)
	            k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        }else if(i != r1 && i != r2){
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        }
	        sample[i].lattice = a_lattice;
	    }
            break;

        case 12:
	    // Place most sequences all in one deme (0,0)
        // and 3 random sequences in random, distinct, demes
        // and Q in (2,2)
	    // d x d = 9
	    // gps_3other_large_rand
            if(debug){ printf("gps 3other large rand\n"); }
	    r1 = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    do{ r2 = SEQ_PER_SPECIES*ran1(iseed); // choose 2nd random seq that has not been chosen
	    }while(r1==r2); // ensure r1 != r2
        do{ r3 = SEQ_PER_SPECIES*ran1(iseed);
        }while(r1==r3 || r2==r3); // r1 != r3 and r2 != r3
	    // do not start from 0 (query)
	    if(begin_nn == 0){
	        r1 = begin_nn + 1 + r1;
	        r2 = begin_nn + 1 + r2;
            r3 = begin_nn + 1 + r3;
	    }else{
	        r1 = begin_nn + r1;
	        r2 = begin_nn + r2;
            r3 = begin_nn + r3;
	    }
        // choose random deme for each random seq
        for(i=0; i<=2; i++){
            do{ // choose random deme
                j = d*ran1(iseed); //choose x
                l = d*ran1(iseed); //choose y
                m = d*j + l; // get deme
            }while((j == 0 && l == 0) || (j == 2 && l == 2) || k_lattice[a_lattice][m] != 0);
            // choose new random deme if is a) base deme, b) deme with query and c) is occupied
            if(i==0){
                sample[r1].deme_x = j;
                sample[r1].deme_y = l;
            }else if (i==1){
                sample[r2].deme_x = j;
                sample[r2].deme_y = l;
            }else if (i==2){
                sample[r3].deme_x = j;
                sample[r3].deme_y = l;
            }
	        k_lattice[a_lattice][m]++;
        } // end random seq
        // place non-random seqs
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = sample[i].deme_y = 2; //(2,2)
	            k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        }else if(i != r1 && i != r2 && i != r3){
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        }
	        sample[i].lattice = a_lattice;
	    }
            break;

        case 13:
	    // Place most sequences all in one deme (0,0)
        // and 4 random sequences in random, distinct, demes
        // and Q in (2,2)
	    // d x d = 9
	    // gps_4other_large_rand
            if(debug){ printf("gps 4other large rand\n"); }
	    r1 = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    do{ r2 = SEQ_PER_SPECIES*ran1(iseed); // choose 2nd random seq that has not been chosen
	    }while(r1==r2); // ensure r1 != r2
        do{ r3 = SEQ_PER_SPECIES*ran1(iseed);
        }while(r1==r3 || r2==r3); // r1 != r3 and r2 != r3
        do{ r4 = SEQ_PER_SPECIES*ran1(iseed);
        }while(r1==r4 || r2==r4 || r3==r4); // r1 != r3 and r2 != r3 and r3 != r4
	    // do not start from 0 (query)
	    if(begin_nn == 0){
	        r1 = begin_nn + 1 + r1;
	        r2 = begin_nn + 1 + r2;
            r3 = begin_nn + 1 + r3;
            r4 = begin_nn + 1 + r4;
	    }else{
	        r1 = begin_nn + r1;
	        r2 = begin_nn + r2;
            r3 = begin_nn + r3;
            r4 = begin_nn + r4;
	    }
        // choose random deme for each random seq
        for(i=0; i<=3; i++){
            do{ // choose random deme
                j = d*ran1(iseed); //choose x
                l = d*ran1(iseed); //choose y
                m = d*j + l; // get deme
            }while((j == 0 && l == 0) || (j == 2 && l == 2) || k_lattice[a_lattice][m] != 0);
            // choose new random deme if is a) base deme, b) deme with query and c) is occupied
            if(i==0){
                sample[r1].deme_x = j;
                sample[r1].deme_y = l;
            }else if (i==1){
                sample[r2].deme_x = j;
                sample[r2].deme_y = l;
            }else if (i==2){
                sample[r3].deme_x = j;
                sample[r3].deme_y = l;
            }else if (i==3){
                sample[r4].deme_x = j;
                sample[r4].deme_y = l;
            }
	        k_lattice[a_lattice][m]++;
        } // end random seq
        // place non-random seqs
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = sample[i].deme_y = 2; //(2,2)
	            k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        }else if(i != r1 && i != r2 && i != r3 && i != r4){
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        }
	        sample[i].lattice = a_lattice;
	    }
            break;

        case 14:
	    // Place most sequences all in one deme (0,0)
            // and 5 random sequences in random, distinct, demes
            // and Q in (2,2)
	    // d x d = 9
	    // gps_5other_large_rand
            if(debug){ printf("gps 5other large rand\n"); }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = sample[i].deme_y = 2; //(2,2)
	            k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        }else{
                do{ // choose unoccupied random deme
                    j = d*ran1(iseed); //choose x
                    l = d*ran1(iseed); //choose y
                    m = d*j + l; // get deme
                }while((j == 2 && l == 2) || k_lattice[a_lattice][m] != 0);
                sample[i].deme_x = j;
                sample[i].deme_y = l;
	            k_lattice[a_lattice][m]++;
	        }
	        sample[i].lattice = a_lattice;
	    }
            break;

        case 15:
	    // Place all sequences in one deme (0,0)
            // and Q in (2,2)
	    // d x d = 9
	    // gps_all_large_rand
            if(debug){ printf("gps all large rand\n"); }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = 2;
	            sample[i].deme_y = 2;
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            sample[i].lattice = a_lattice;
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        case 16:
	    // Place all sequences in one deme (0,0)
            // and Q in (3,3)
	    // d x d = 16
	    // gps_all_large_4x4
            if(debug){ printf("gps all large 4x4\n"); }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = 3;
	            sample[i].deme_y = 3;
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            sample[i].lattice = a_lattice;
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;
	    
        case 17:
	    // Place most sequences all in one deme (0,0)
	    // except for 1 randomly chosen sequence (3,3) and
	    // Q in (3,0)
	    // d x d = 16
	    // gps_1other_large_4x4
            if(debug){ printf("gps 1other large 4x4\n"); }
	    j = SEQ_PER_SPECIES*ran1(iseed); // choose a random seq on a scale of 0-4
	    // the chosen random seq
	    if(begin_nn == 0){ j = begin_nn + 1 + j; } // do not start from 0 (query)
	    else{ j = begin_nn + j; }
	    for(i=begin_nn; i<=end_nn; i++){
	        if(i==0){
	            sample[i].deme_x = 3;
	            sample[i].deme_y = 0;
	        }else if(i==j){
	            sample[i].deme_x = sample[i].deme_y = 3; //(3,3)
	        }else{
	            sample[i].deme_x = sample[i].deme_y = 0; //(0,0)
	            sample[i].lattice = a_lattice;
	        }
	        sample[i].lattice = a_lattice;
	        k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	    }
            break;

        default:
	    // Place sequences randomly in a k_lattice[a_lattice][d*d]
            if(debug){ printf("default: set randomly\n"); }
	    for(i=begin_nn; i<=end_nn; i++){
	        j = d*ran1(iseed);
	        sample[i].deme_x = j;
		j = d*ran1(iseed);
	        sample[i].deme_y = j;
		sample[i].lattice = a_lattice; // place in new lattice
		k_lattice[a_lattice][d*sample[i].deme_x + sample[i].deme_y]++;
	        if(debug){ printf("%d ",i); }
	    }
            break;
    }//switch
}
/********************************************************************/
void   set_uncoalesced_position(int a_seq, int a_lattice, int d, double scale) {
    int i, j;
    int debug = 0;
    int scheme = 0; // default
    /* data structure:
    int    deme_x, deme_y;      // deme (x,y) coordinates in demes number
    double gps_x, gps_y;        // deme (w,z) coordinates in GPS hard values
    */
    //-----------------------------------------------------------------------
    switch(scheme){
        case 1:
        // all: place uncoalesced sequence(s) in (0,0)
            if(sample[a_seq].anc == 0 && sample[a_seq].lattice != a_lattice){
                if(debug){ printf("\nall: set uncoalesced ... %d in lattice %d\n", a_seq, a_lattice); }
                sample[a_seq].deme_x = sample[a_seq].deme_y = 0; // (0,0)
                sample[a_seq].lattice = a_lattice; // place in new lattice
                k_lattice[a_lattice][d*sample[a_seq].deme_x + sample[a_seq].deme_y]++;
            }
        break;
        // random:
        // place uncoalesced sequence (anc == 0 && anc_time == begin_time (see
        // do_block_coalescent) && lattice != a_lattice) randomly in a
        // k_lattice[a_lattice][d*d]
        default:
            if(sample[a_seq].anc == 0 && sample[a_seq].lattice != a_lattice){
                if(debug){ printf("\nrandom: set uncoalesced ... %d in lattice %d\n", a_seq, a_lattice); }
                j = d*ran1(iseed);
                sample[a_seq].deme_x = j;
                j = d*ran1(iseed);
                sample[a_seq].deme_y = j;
                sample[a_seq].lattice = a_lattice; // place in new lattice
                k_lattice[a_lattice][d*sample[a_seq].deme_x + sample[a_seq].deme_y]++;
            }
        break;
    }// end switch
}
/********************************************************************/
// choose a gene from deme d from total of k sequences given the total
// number demes, demes, and sequences, nn, on the lattice
int choose_a_gene(int d, int k, int demes, int nn, int a_lattice) {
    int i, j, count=-1, current_deme;
    int debug = 0;
    int seq;
    
    j = ran1(iseed)*k; // choose a random sequence of k sequences from deme d
    if (debug){
        printf("from pop d=%d, picking the %dth seq out of %d seqs in %d demes ",d, j, k, demes);
        printf("from lattice %d\n", a_lattice);
        for(i=0; i<nn; i++){
            seq = k_seq[a_lattice][i];
	    printf("Seq %d:\t", seq);
	    printf("(d_x, d_y): (%d, %d)\t", sample[seq].deme_x, sample[seq].deme_y);
            printf("deme: %d\n", demes*sample[seq].deme_x + sample[seq].deme_y);
        }
    }
    for(i=0; i<nn; i++){
        seq = k_seq[a_lattice][i];
	current_deme = demes * sample[seq].deme_x + sample[seq].deme_y;
	if(current_deme == d) count++; // in deme, d
        if (debug){ printf("seq: %d curr deme: %d cnt: %d\n", seq, current_deme, count); }
	if(j==count) return(seq);
    }
    fprintf(stderr,"ERROR in choose_a_gene, from pop d=%d, picking the %dth seq ",d,j); // shouldn't reach this line
    fprintf(stderr," from a total of %d seqs, demes=%d total  \n",k,demes); 
    fprintf(stderr," from lattice %d\n", a_lattice);
    exit(0);
}
/********************************************************************/
// Find descendants of node i 
// 2n-1 equals the exact number of nodes required to build a tree
a_desc  find_descendants(int a_seq) { 
    int a,b,j;
    a_desc desc;
    a=b=-1;
    for(j=0; j<number_of_seq*2 - 1; j++) {
	if(sample[j].anc == a_seq) break;
    }
    a=j; // descendant 1 of a_seq
    for(j=a+1; j<number_of_seq*2 - 1; j++) {
	if(sample[j].anc == a_seq) break;
    }
    b=j; // descendant 2 of a_seq
    if(a==-1 || b==-1) { 
	printf(" Error finding descendants of i=%d\n",a_seq);
	exit(0);
    }
    desc.desc1=a;
    desc.desc2=b;
    return desc;
}
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/

/* Get in-order tree */
void inorder(int n){
    // n holds the root of the tree and represents the number of nodes
    // in the array
    // n+1 holds the actual number of nodes
    int debug = 0;

    trackAncCount = 0;
    newickCount = 0;
    memset(newickTree, '\0', sizeof(newickTree));
    memset(trackAnc, 0, sizeof(trackAnc)); // initialize array to track visited ancestor nodes
    newickRoot = n;
    inorder_left(n);
    
    if (debug){ printf("%s\n", newickTree); }// print final newick tree
    //printf("%d\n", newickCount); // count chars used for tree
}

/* Print left most node */
void inorder_left(int n){
    char buffer[30];
    // assume the taxa_id and its distance is not represented by more than 30 characters
    int i = 0;
    int debug = 0;
    int leftNode;
    int leftleftNode;
    int j = 0;
    
    leftNode = nj_tree[n].left;
    leftleftNode = nj_tree[leftNode].left;

    //debug
    if (debug){
        printf("current node: %d\n", nj_tree[n].taxa_id);
        printf("current node's left node: %d\n", leftNode);
        printf("current node's left node left node: %d\n", leftleftNode);
        for (j = 0; j <= newickRoot; j++) printf("LN trackAnc[%d] = %d\n", j, trackAnc[j]);
    }

    while (leftleftNode != -1){ //left left node
        trackAnc[trackAncCount++] = n;
        n = nj_tree[n].left;
        leftNode = nj_tree[n].left;
        leftleftNode = nj_tree[leftNode].left;
        if (debug){
            printf("LN1\n");
            printf("current node: %d\n", nj_tree[n].taxa_id);
            printf("current node's left node's: %d\n", leftNode);
            printf("current node's left node left node: %d\n", leftleftNode);
        }
        newickTree[newickCount++] = '(';
    }
    trackAnc[trackAncCount++] = n;
    n = nj_tree[n].left;
    
    // debug print ancestor array
    if (debug){ for (j = 0; j <= newickRoot ; j++) printf("trackAnc[%d] = %d\n", j, trackAnc[j]); }
    
    // record the leftest node
    newickTree[newickCount++] = '(';
    // convert integer to string
    memset(buffer, '\0', 30); // clear the array
    //sprintf(buffer, "%d", nj_tree[n].taxa_id);
    sprintf(buffer, "%s", sample[nj_tree[n].taxa_id].title);
    //printf("title of nj_tree[%d] with taxa id %d: %s\n", n, nj_tree[n].taxa_id, sample[nj_tree[n].taxa_id].title);
    // append name
    // append each digit (now a char) in buffer to newickTree
    i = 0;
    while (buffer[i] != '\0'){ newickTree[newickCount++] = buffer[i++]; }
    if (debug){ printf("\n"); }
    // append distance
    newickTree[newickCount++] = ':';
    memset(buffer, '\0', 30);
    sprintf(buffer, "%f", nj_tree[n].distance);
    
    i=0;
    while (buffer[i] != '\0'){ newickTree[newickCount++] = buffer[i++]; }

    if (debug){
        printf("LN2 newick tree: %s\n", newickTree);
    }
    
    inorder_right(trackAnc[trackAncCount-1]); // go back up to ancestor of leftest node
}

/* Print right most node */
void inorder_right(int n){
    char buffer[30];
    int i = 0;
    int debug = 0;
    int j= 0;
    int rightNode;

    rightNode = nj_tree[n].right;
    if (debug){
        printf("current node: %d\n", nj_tree[n].taxa_id);
        printf("current node's right node: %d\n", nj_tree[n].right);
    }
    if (nj_tree[n].checked != 1){
        nj_tree[n].checked = 1; // mark node checked when all left node(s) visited
        if (debug){ printf("checked node %d\n", n); }
        if (nj_tree[rightNode].left != -1){ // if the right node has a left node
            if (debug){ printf("RN1\n"); }
            newickTree[newickCount++] = ',';
            inorder_left(rightNode); // use the right node
        }else{ // found the right most node of current ancestor
            newickTree[newickCount++] = ',';
            memset(buffer, '\0', 30); // clear the array
            // append name
            //sprintf(buffer, "%d", nj_tree[rightNode].taxa_id);
            sprintf(buffer, "%s", sample[nj_tree[rightNode].taxa_id].title);
            i = 0;
            while (buffer[i] != '\0'){ newickTree[newickCount++] = buffer[i++]; }
            // append distance
            newickTree[newickCount++] = ':';
            memset(buffer, '\0', 30);
            sprintf(buffer, "%f", nj_tree[rightNode].distance);
            i = 0;
            while (buffer[i] != '\0'){ newickTree[newickCount++] = buffer[i++]; }
            if (debug){ printf("RN2 newick tree: %s\n", newickTree); }
            inorder_up(n);
        }//end else
    }
}

/* Print ancestor node whose decendents have been printed */
void inorder_up(int n){
    char buffer[30];
    int i = 0;
    int debug = 0;
    int j= 0;

    while (nj_tree[n].checked == 1){
        if (nj_tree[n].taxa_id == newickRoot){
            if (debug){ printf("FINAL!\n"); }
            newickTree[newickCount++] = ')';
            newickTree[newickCount++] = ';';
        }else{
            newickTree[newickCount++] = ')';
            newickTree[newickCount++] = ':';
            memset(buffer, '\0', 30);
            sprintf(buffer, "%f", nj_tree[n].distance);
            i = 0;
            while (buffer[i] != '\0'){ newickTree[newickCount++] = buffer[i++]; }
            if (debug){ printf("UN2 newick tree: %s\n", newickTree); }
            // remove current node from trackAnc
            if(debug){ for (j = 0; j <= newickRoot; j++) printf("UN2 trackAnc[%d] = %d\n", j, trackAnc[j]); }
            trackAnc[--trackAncCount] = '\0';
            if(debug){ for (j = 0; j <= newickRoot; j++) printf("UN2 trackAnc[%d] = %d\n", j, trackAnc[j]); }
        }
        if (debug){ printf("UN POST trackAnc[%d] = %d\n", trackAncCount-1, trackAnc[trackAncCount-1]); }
        if (n == newickRoot && nj_tree[n].checked == 1){ break;
        }else{
            n = trackAnc[trackAncCount-1];
            if(debug) {printf("if not checked-root, is %d checked?\n", n); }
        }
    }
    inorder_right(n);
}
/********************************************************************/
void print_sample(int i, int d){
    printf("id: %d\tl: %d\t", i, sample[i].lattice);
    printf("(%d, %d)\t", sample[i].deme_x, sample[i].deme_y);
    printf("%5d\t", d*sample[i].deme_x+sample[i].deme_y);
    printf("\ta: %d\t(a_t): %.5f\t%.5f\t", sample[i].anc, sample[i].anc_time, sample[i].time);
    printf("\td1: %d\td2: %d\n", sample[i].desc1, sample[i].desc2);
    printf("mutns: %d\n",sample[i].mutations);
    printf("visited: l=%d\tr=%d\n",sample[i].l_visited,sample[i].r_visited);
}
/********************************************************************/
void print_block(int begin_nn, int end_nn, int a_lattice, int d, int nn){
    int i;

    if(print2File){
        fprintf(out_debug,"print block ... %d: %d - %d\n", a_lattice, begin_nn, end_nn);
    }else{
        printf("print block ... %d: %d - %d\n", a_lattice, begin_nn, end_nn);
    }
    // print k_seq (seqs not undergone an event)
    for(i=0; i<nn; i++){
        if(print2File){ fprintf(out_debug,"     uncoal k_seq[%d][%d]: %d\n", a_lattice, i, k_seq[a_lattice][i]); }
        else{ printf("     uncoal k_seq[%d][%d]: %d\n", a_lattice, i, k_seq[a_lattice][i]); }
    }
    // print info for each seq in variable, sample
    for(i=begin_nn; i<=end_nn; i++){
        if(sample[i].lattice == a_lattice){
            if(print2File){
                fprintf(out_debug,"%6d:%5d   ", i, sample[i].lattice);
                fprintf(out_debug,"(%3d,%3d)", sample[i].deme_x, sample[i].deme_y);
                fprintf(out_debug,"%5d\t", d*sample[i].deme_x+sample[i].deme_y);
                fprintf(out_debug,"a: %d\t(a_t): %.5f\tn_t: %.5f\t", sample[i].anc, sample[i].anc_time, sample[i].time);
                fprintf(out_debug,"d1: %d\td2: %d\n", sample[i].desc1, sample[i].desc2);
            }else{
               printf("%6d:%5d   ", i, sample[i].lattice);
               printf("(%3d,%3d)", sample[i].deme_x, sample[i].deme_y);
               printf("%5d\t", d*sample[i].deme_x+sample[i].deme_y);
               printf("a: %d\t(a_t): %.5f\tn_t: %.5f\t", sample[i].anc, sample[i].anc_time, sample[i].time);
               printf("d1: %d\td2: %d\n", sample[i].desc1, sample[i].desc2);
            }
        }
    }
    for(i=MAX_SEQ; i<old_anc; i++) {
        if(sample[i].lattice == a_lattice){
            if(print2File){
                fprintf(out_debug,"%6d:%5d   ", i, sample[i].lattice);
                fprintf(out_debug,"(%3d,%3d)", sample[i].deme_x, sample[i].deme_y);
                fprintf(out_debug,"%5d\t", d*sample[i].deme_x+sample[i].deme_y);
                fprintf(out_debug,"a: %d\t(a_t): %.5f\tn_t: %.5f\t", sample[i].anc, sample[i].anc_time, sample[i].time);
                fprintf(out_debug,"d1: %d\td2: %d\n", sample[i].desc1, sample[i].desc2);
            }else{
               printf("%6d:%5d   ", i, sample[i].lattice);
               printf("(%3d,%3d)", sample[i].deme_x, sample[i].deme_y);
               printf("%5d\t", d*sample[i].deme_x+sample[i].deme_y);
               printf("a: %d\t(a_t): %.5f\tn_t: %.5f\t", sample[i].anc, sample[i].anc_time, sample[i].time);
               printf("d1: %d\td2: %d\n", sample[i].desc1, sample[i].desc2);
            }
        }
    }
    for(i=old_anc; i<max_anc; i++) {
        if(sample[i].lattice == a_lattice){
            if(print2File){
                fprintf(out_debug,"%6d:%5d   ", i, sample[i].lattice);
                fprintf(out_debug,"(%3d,%3d)", sample[i].deme_x, sample[i].deme_y);
                fprintf(out_debug,"%5d\t", d*sample[i].deme_x+sample[i].deme_y);
                fprintf(out_debug,"a: %d\t(a_t): %.5f\tn_t: %.5f\t", sample[i].anc, sample[i].anc_time, sample[i].time);
                fprintf(out_debug,"d1: %d\td2: %d\n", sample[i].desc1, sample[i].desc2);
            }else{
               printf("%6d:%5d   ", i, sample[i].lattice);
               printf("(%3d,%3d)", sample[i].deme_x, sample[i].deme_y);
               printf("%5d\t", d*sample[i].deme_x+sample[i].deme_y);
               printf("a: %d\t(a_t): %.5f\tn_t: %.5f\t", sample[i].anc, sample[i].anc_time, sample[i].time);
               printf("d1: %d\td2: %d\n", sample[i].desc1, sample[i].desc2);
            }
        }
    }
    if(print2File){ fprintf(out_debug,"finish block\n"); }
    else{ printf("finish block\n"); }
}
/********************************************************************/
void print_lattice(int a_lattice, int d){
    int i,sum;

    printf("print seqs and locations for lattice: %d ...\n", a_lattice);
    for(i=0,sum=0; i<d*d; i++){
        printf("%d ", k_lattice[a_lattice][i]);
        sum+=k_lattice[a_lattice][i];
    }
    printf(": %d\n",sum);
    for(i=0; i<MAX_NODES; i++){
        if (sample[i].lattice == a_lattice){
            printf("\tseq id: %d deme: %d\n", i, d*sample[i].deme_x + sample[i].deme_y);
        }
    }
}
/********************************************************************/
/* find the largest seq (from a different taxon) that coalesces with 0
 */
int find_first(int a_seq) { // find first coalescent with a_seq
    int i, j, m, k[MAX_SEQ], count=0;
    int anc;
    int debug = 0;
    
    anc = sample[a_seq].anc;
    // record descendants of a_seq
    k[count++]=sample[anc].desc1; 
    k[count++]=sample[anc].desc2; 
    if(debug){
        if(print2File){
            fprintf(out_debug,"find_first(%d)\n",a_seq);
            fprintf(out_debug,"\tanc(%d): ", anc);
            fprintf(out_debug,"%d %d\n", sample[anc].desc1, sample[anc].desc2);
        }else{
            printf("find_first(%d)\n",a_seq);
            printf("\tanc(%d): ", anc);
            printf("%d %d\n", sample[anc].desc1, sample[anc].desc2);
        }
    }
    i=0;
    // if it's an internal node (i.e. > MAX_SEQ), then record its'
    // descendants too
    do {
	if(k[i]>=MAX_SEQ) {
	    j=k[i];
	    k[i]=sample[j].desc1;  // replace original
	    k[count++]=sample[j].desc2; 
	    i=0;
	}
	else { i++; }
    } while(i < count);
    // find largest k[]; found at pos 0
    // k
    for(j=1; j<count; j++) {
	if(k[0]<k[j]) {
	    m=k[j];
	    k[j]=k[0];
	    k[0]=m;
	}
    }
    if(debug){
        if(print2File){
            fprintf(out_debug,"\t");
            for(i=0; i<count; i++){ fprintf(out_debug,"%d ", k[i]); }
            fprintf(out_debug,"\n");
        }else{
            printf("\t");
            for(i=0; i<count; i++){ printf("%d ", k[i]); }
            printf("\n");
        }
    }
    return(k[0]); // return largest value less than MAX_SEQ
}
/********************************************************************/
void find_all(int a_seq) { // find all descendants of node a_seq
    int i, j, count=0;
    int debug = 0;

    for(i=0; i<MAX_SEQ; i++){ k[i] = -1; } // clear
    k[count++]=sample[a_seq].desc1; 
    k[count++]=sample[a_seq].desc2; 
    i=0;
    do {
	if(k[i]>=MAX_SEQ) {
	    j=k[i];
	    k[i]=sample[j].desc1;  // replace original with desc1
	    k[count++]=sample[j].desc2; 
	    i=0;
	}
	else { i++; }
    } while(i < count);
    if(debug){
        if(print2File){
            fprintf(out_debug,"\tfind_all(%d): ",a_seq);
            for(i=0; i<count; i++){ fprintf(out_debug,"%d ", k[i]); }
            fprintf(out_debug,"\n");
        }else{
            printf("\tfind_all(%d): ",a_seq);
            for(i=0; i<count; i++){ printf("%d ", k[i]); }
            printf("\n");
        }
    }
}
/********************************************************************/
int find_all_taxa1(int n) {  
    // find last common ancestor of sequences from first taxa (taxa0)
    // excluding seq0 (the query) and return the no of sequences from
    // other taxa that are included in the coalescent with seqs 1..5
    int i, j, m, found[SEQ_PER_SPECIES+1];
    int debug = 0;
    
    if(debug){
        if(print2File){ fprintf(out_debug,"\nfind_all_taxa1(%d)\n", n); }
        else{ printf("\nfind_all_taxa1(%d)\n", n); }
    }
    n = sample[n].anc;
    if(debug){
        if(print2File){ fprintf(out_debug,"\tanc: %d\n", n); }
        else{ printf("\tanc: %d\n", n); }
    }
    find_all(n); // find all descendants of node n
    // are 1-5 all in k?
    for(i=0; i<SEQ_PER_SPECIES+1; i++) {
	found[i] = 0; // initialize
    }
    m = 0;
    if(debug){
        if(print2File){ fprintf(out_debug,"found: ");
        }else{ printf("found: "); }
    }
    for(j=1; j<SEQ_PER_SPECIES+1; j++) { // for each 1-5
	for(i=0; i<MAX_SEQ; i++) { // for each 0 .. 50
            if(k[i] == j+m){
                found[j] = 1;
                if(debug){
                    if(print2File){ fprintf(out_debug,"%d ",j); }
                    else{ printf("%d ",j); }
                }
            }
        }
    }
    if(debug){
        if(print2File){ fprintf(out_debug,"\n"); }
        else{ printf("\n"); }
    }
    // yes 0-5 are all in k ... return n
    // needs to actually return largest subtended by n
    for(j=1, i=1; i<SEQ_PER_SPECIES+1; i++) {
	j=j*found[i]; 
    }
    if(j == 1) return(n);
    // no  1-5 are not in k ... go around
    n = find_all_taxa1(n);
}
/********************************************************************/
void print_mutations(int a_seq){
    if(a_seq == MAX_SEQ*2-2){
        if(print2File){ fprintf(out_debug,"%5d: %d\n",a_seq,sample[a_seq].mutations); }
        else{ printf("%5d: %d\n",a_seq,sample[a_seq].mutations); }
    }else{
        if(print2File){ fprintf(out_debug,"%5d: %d",a_seq,sample[a_seq].mutations); }
        else{ printf("%5d: %d",a_seq,sample[a_seq].mutations); }
        print_mutations(sample[a_seq].anc);
    }
}
/********************************************************************/
int find_last(int a_taxon, int a_seq) {
    // find common ancestor of seqs in a_taxon (starting with a_seq)
    // list all descendants of seq a_seq, then find lowest number that
    // includes all of seqs from a_seq ... a_seq+SEQ_PER_SEQ
    int i, found1, found2, found3, found4, found5;
    int debug = 0;
    int startSeq = a_taxon*SEQ_PER_SPECIES+1;
    
    if(debug){
        if(print2File){ fprintf(out_debug,"find_last(%d)\n", a_seq); }
        else{ fprintf(stdout, "find_last(%d)\n", a_seq); }
    }
    a_seq = sample[a_seq].anc; // get its' ancestor
    if(debug){
        if(print2File){ fprintf(out_debug,"\tanc: %d\n", a_seq); }
        else{ fprintf(stdout, "find_last(%d)\n", a_seq); }
    }
    // find all descendants of a_seqs' ancestor and place in global
    // variable k[]
    find_all(a_seq);
    // are 0-5 all in k?
    found1 = 0;
    found2 = 0;
    found3 = 0;
    found4 = 0;
    found5 = 0;
    for(i=0; i<MAX_SEQ; i++) {
	if(k[i] == startSeq) found1 = 1;
	if(k[i] == startSeq+1) found2 = 1;
	if(k[i] == startSeq+2) found3 = 1;
	if(k[i] == startSeq+3) found4 = 1;
	if(k[i] == startSeq+4) found5 = 1;
    }
    // yes 0-5 are all in k ... return a_seq
    // needs to actually return largest subtended by a_seq
    if(found1*found2*found3*found4*found5 == 1){ return(a_seq); }
    // no  0-5 are not in k ... check its' ancestor (recursive)
    a_seq = find_last(a_taxon,a_seq);
}
/********************************************************************/
int find_all_incl(int a_seq) { // find all descendants of node a_seq
    int i, j, count=0;
    int debug = 0;
    
    for(i=0; i<MAX_SEQ; i++){ k[i] = -1; } // clear
    k[count++]=sample[a_seq].desc1; 
    k[count++]=sample[a_seq].desc2; 
    i=0;
    do {
	if(k[i]>=MAX_SEQ) {
	    j=k[i];
	    k[count++]=sample[j].desc1;
	    k[count++]=sample[j].desc2; 
	}
        i++;
    } while(i < count);
    
    if(debug){
        if(print2File){
            fprintf(out_debug,"find_all_incl(%d): ",a_seq);
            for(i=0; i<count; i++){ fprintf(out_debug,"%d ", k[i]); }
            fprintf(out_debug,"\n");
        }else{
            printf("find_all_incl(%d): ",a_seq);
            for(i=0; i<count; i++){ printf("%d ", k[i]); }
            printf("\n");
        }
    }
    return count;
}
/********************************************************************/
/* Transform dxd to k */
void print_transform(int d) {
    int i,j,l;

    for(i=0; i<d; i++){
        for(j=0;j<d;j++){
            l=(d*i)+j;
            printf("(%d, %d): %d\n",i,j,l);
        }
    }
}

