#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <regex>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <random> // shuffle
#include <chrono> // chrono
#include <algorithm> // upper_bound, sort
#include <regex.h>

using namespace std;

#define MAX_ERROR_MSG 0x1000

/* Global Variables */
static int permute_round;
static vector< vector< unordered_set<int> > > deg_list; // 3D deg_list to store DEG lists for each spe, 0:hb, 1:mm, 2:sb [spe1][deg_list1][gID1], where gID1 is the pos in arr
//static vector< vector< unordered_set<int> > > deg_list_int; // 3D deg_list to store DEG lists for each spe, 0:hb, 1:mm, 2:sb [spe1][deg_list1][gID1], this one should be obtained after read in orth file
static vector< vector<string> > deg_list_name (3, vector<string>(0)); // 2D arr to store deg_list file names
static vector< vector<string> > universe_gene_set(3, vector<string>(0)); // 2D arr to store unique ortholog genes for each spes: [spe1][gID1], this arr can be permuted (swap gID) to generate random permutation for enrich test
static unordered_set<string> orth_edges; // 4D arr to store ortholog edges between each pair of spes: [spe1][spe2][gID1][gID2]
static unordered_map<string, int> gNameToID; // mapping from geneName to geneID
//static vector< unordered_map<int, string> > gIDToName; // mapping from geneID to geneName
static vector<string> deg_pair; // arr to store all deg pairs
static vector< vector<int> > enrich_deg; // 2D arr to store permut orth edge count of deg pairs
static vector< int > enrich_deg_origUniverse; // arr to store orig orth edge count of deg pairs
static vector< int > deg_rank; // arr to store the rank of deg in permut distribution

/* Function Declaration */
void readDEGfile(char* input_dir, char* deg_list_file);
void readOrthGroup(char* orth_file);
int orthCount(unordered_set<int>& deg_set1, unordered_set<int>& deg_set2, vector<string>& universe1, vector<string>& universe2);
vector<string> permutUniverse(vector<string>& orig_arr);
void enrichInPermutUniverse();
void readUniverse(char* universe_file);
int compile_regex (regex_t * r, char * regex_text);
int match_regex (regex_t * r, char * to_match);

int compile_regex (regex_t * r, char * regex_text)
{
    int status = regcomp (r, regex_text, REG_EXTENDED|REG_NEWLINE);
    if (status != 0) {
        char error_message[MAX_ERROR_MSG];
        regerror (status, r, error_message, MAX_ERROR_MSG);
        printf ("Regex error compiling '%s': %s\n",
                regex_text, error_message);
        return 1;
    }
    return 0;
}

int match_regex (regex_t * r, char * to_match)
{
    /* "P" is a pointer into the string which points to the end of the
     previous match. */
    const char * p = to_match;
    /* "N_matches" is the maximum number of matches allowed. */
    const int n_matches = 10;
    /* "M" contains the matches found. */
    regmatch_t m[n_matches];
    
    int nomatch = regexec (r, p, n_matches, m, 0);
    if (!nomatch)
    {
        return 1;
    }
    return 0;
}


void rankOfDGE()
{
    for (int i = 0; i<enrich_deg.size();i++)
    {
        sort(enrich_deg[i].begin(), enrich_deg[i].end());
        fprintf(stderr, "deg pair\t%d: ", i);
        for (auto x:enrich_deg[i]) fprintf(stderr, "%d ", x);
        fprintf(stderr, "\n");
        vector<int>::iterator low;
        low = lower_bound(enrich_deg[i].begin(), enrich_deg[i].end(), enrich_deg_origUniverse[i]);
        deg_rank.push_back((int)(enrich_deg[i].end() - low));
    }
    
    for (int i = 0; i<deg_rank.size();i++)
    {
        printf("rank of\tdeg_pair\t%s\t%d\n", deg_pair[i].c_str(), deg_rank[i]);
    }
}


int orthCount(unordered_set<int>& deg_set1, unordered_set<int>& deg_set2, vector<string>& universe1, vector<string>& universe2)
{
    int orth_edge_count = 0;
    unordered_set<int>::iterator itr1, itr2;
    for (itr1 = deg_set1.begin(); itr1!=deg_set1.end();itr1++)
    {
        if (*itr1 < 0 || *itr1 > universe1.size())
        {
            printf("invalid pos %d, with universe size %lu\n", *itr1, universe1.size());
            exit(1);
        }
        string g1 = universe1[*itr1];
        
        //        printf("gene1 pos %d, gene name %s\n", *itr1, g1.c_str());
        for (itr2 = deg_set2.begin(); itr2!=deg_set2.end(); itr2++)
        {
            string g2 = universe2[*itr2];
            //            printf("gene2 pos %d, gene name %s\n", *itr2, g2.c_str());
            string gene_pair = g1 + "_" + g2;
            if (orth_edges.find(gene_pair) != orth_edges.end())
            {
                orth_edge_count++;
            }
        }
    }
    return orth_edge_count;
}

vector<string> permutUniverse(vector<string>& orig_arr)
{
    vector<string> permut_arr(orig_arr);
    /* random shuffle arr */
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    auto gen = std::default_random_engine(seed);
    shuffle(permut_arr.begin(), permut_arr.end(), gen);
    return permut_arr;
}

void enrichInPermutUniverse()
{
    /* generate all possible combination of deg pairs */
    int deg_pair_counter = 0;
    for (int i = 0; i<3;i++)
    {
        for (int j = i+1;j<3;j++)
        {
            for (int m = 0; m<deg_list_name[i].size();m++)
            {
                for (int n = 0; n<deg_list_name[j].size(); n++)
                {
                    string tmp_deg_pair = deg_list_name[i][m] + " " + deg_list_name[j][n];
                    deg_pair.push_back(tmp_deg_pair);
                    vector<int> tmp_deg_permut_test;
                    enrich_deg.push_back(tmp_deg_permut_test);
                    /* compute orth count for deg pairs with orig universe */
                    int orth_count = orthCount(deg_list[i][m], deg_list[j][n], universe_gene_set[i], universe_gene_set[j]);
                    fprintf(stderr, "orig orth\t%s\tcount\t%d\n", deg_pair[deg_pair_counter].c_str(), orth_count);
                    enrich_deg_origUniverse.push_back(orth_count);
                    deg_pair_counter++;
                }
            }
        }
    }
    /* do 10000 permut test */
    for (int k = 0; k<permute_round; k++)
    {
        fprintf(stderr, "permut round %d\n", k);
        /* each iter permut all spes universe gene set */
        vector< vector<string> > permut_universe(3, vector<string>(0));
        for (int i = 0; i<3;i++)
        {
            permut_universe[i] = permutUniverse(universe_gene_set[i]);
            // for (auto x:permut_universe[i]) printf("%s ", x.c_str());
            // printf("\n");
            //            printf("size of permut universe %d\n", permut_universe[i].size());
        }
        /* count orth edges for each pair of degs in permut universe genes */
        deg_pair_counter = 0;
        for (int i = 0; i<3;i++)
        {
            for (int j = i+1;j<3;j++)
            {
                for (int m = 0; m<deg_list[i].size();m++)
                {
                    for (int n = 0; n<deg_list[j].size(); n++)
                    {
                        int orth_count = orthCount(deg_list[i][m], deg_list[j][n], permut_universe[i], permut_universe[j]);
                        // fprintf(stderr, "permute orth count %d\n", orth_count);
                        enrich_deg[deg_pair_counter].push_back(orth_count);
                        deg_pair_counter++;
                    }
                }
            }
        }
    }
}

/* read in universe gene set */
void readUniverse(char* universe_file)
{
    int hb_geneID_counter = 0, mm_geneID_counter = 0, sb_geneID_counter = 0;
    string gene;
    smatch m;
    regex_t HB, MM, SB;
    compile_regex(& HB, "(GB)(.*)");
    compile_regex(& MM, "(ENSMUSG)(.*)");
    compile_regex(& SB, "(ENSGACG)(.*)");
    ifstream UNIFILE(universe_file);
    if (UNIFILE.is_open())
    {
        while (getline(UNIFILE, gene))
        {
//            if (regex_search(gene, m, regex("(GB)(.*)")))
            if (match_regex(& HB, &gene[0]))
            {
                universe_gene_set[0].push_back(gene);
                gNameToID[gene] = hb_geneID_counter;
                hb_geneID_counter++;
            }
//            else if (regex_search(gene, m, regex("(ENSMUSG)(.*)")))
            else if (match_regex(& MM, &gene[0]))
            {
                universe_gene_set[1].push_back(gene);
                gNameToID[gene] = mm_geneID_counter;
                mm_geneID_counter++;
            }
//            else if (regex_search(gene, m, regex("(ENSGACG)(.*)")))
            else if (match_regex(& SB, &gene[0]))
            {
                universe_gene_set[2].push_back(gene);
                gNameToID[gene] = sb_geneID_counter;
                sb_geneID_counter++;
            }
        }
    }
    else
    {
        fprintf(stderr, "cannot open universe file\n");
    }
    UNIFILE.close();
}

/* read in orth groups and hash orth gene pairs */
void readOrthGroup(char* orth_file)
{
    string header, line;
    ifstream ORTHFILE(orth_file);
    getline(ORTHFILE, header);
    if (ORTHFILE.is_open())
    {
        while (getline(ORTHFILE, line))
        {
            vector< vector<string> > orth_group(3, vector<string>(0)); // buff arr to store genes in the same orth group
            istringstream split_line (line);
            string groupID, hb_g, mm_g, sb_g;
            split_line >> groupID >> hb_g >> mm_g >> sb_g;
            string gene;
            istringstream split_hb (hb_g);
            while (getline(split_hb, gene, ';'))
            {
                orth_group[0].push_back(gene);
            }
            istringstream split_mm (mm_g);
            while (getline(split_mm, gene, ';'))
            {
                orth_group[1].push_back(gene);
            }
            istringstream split_sb (sb_g);
            while (getline(split_sb, gene, ';'))
            {
                orth_group[2].push_back(gene);
            }
            /* pair up orth genes */
            for (int i = 0; i<3;i++)
            {
                for (int j = i+1; j<3; j++)
                {
                    // pair g1_g2 as key for set
                    for (int m = 0; m < orth_group[i].size();m++)
                    {
                        for (int n = 0; n < orth_group[j].size();n++)
                        {
                            // add g1_g2 pair to orth_edges set
                            string gene_pair = orth_group[i][m] + "_" + orth_group[j][n];
                            orth_edges.emplace(gene_pair);
                        }
                    }
                }
            }
        }
    }
    else
    {
        fprintf(stderr, "cannot open orth group file\n");
    }
    ORTHFILE.close();
}

/* read in DEG lists from all spes */
void readDEGfile(char* input_dir, char* deg_list_file)
{
    // init deg_list
    for(int i = 0; i<3; i++)
    {
        vector< unordered_set<int> > tmp_spe_deg;
        deg_list.push_back(tmp_spe_deg);
    }
    /* read deg list file */
    string INPUTDIR(input_dir);
    string line;
    ifstream DEGLIST (deg_list_file);
    smatch m;
    regex_t hb, mm, sb;
    compile_regex(& hb, "(hb_)(.*)");
    compile_regex(& mm, "(mm_)(.*)");
    compile_regex(& sb, "(sb_)(.*)");
    if (DEGLIST.is_open())
    {
        while (getline(DEGLIST, line))
        {
            //            printf("%s\n", line.c_str());
            unordered_set<int> tmp_deg;
            
            /* push deg list into deg_list and deg_list_name in same order */
//            if (regex_search(line, m, regex("(hb_)(.*)")))
            if (match_regex(& hb, &line[0]))
            {
                // [0].push_back(deg set)
                deg_list[0].push_back(tmp_deg);
                // [0].push_back(deg file name)
                deg_list_name[0].push_back(line);
            }
            else if (match_regex(& mm, &line[0]))
//            else if (regex_search(line, m, regex("(mm_)(.*)")))
            {
                deg_list[1].push_back(tmp_deg);
                deg_list_name[1].push_back(line);
            }
//            else if (regex_search(line, m, regex("(sb_)(.*)")))
            else if (match_regex(& sb, &line[0]))
            {
                deg_list[2].push_back(tmp_deg);
                deg_list_name[2].push_back(line);
            }
            /* read deg file */
            ifstream DEGFILE (INPUTDIR + "/" + line);
            string deg_gene;
            int deg_count = 0;
            regex_t HB, MM, SB;
            compile_regex(& HB, "(GB)(.*)");
            compile_regex(& MM, "(ENSMUSG)(.*)");
            compile_regex(& SB, "(ENSGACG)(.*)");
            if (DEGFILE.is_open())
            {
                while(getline(DEGFILE, deg_gene))
                {
                    deg_count++;
//                    if (regex_search(deg_gene, m, regex("(GB)(.*)")))
                    if (match_regex(& HB, &deg_gene[0]))
                    {
                        // [0][deg_listID].emplace(gID)
                        deg_list[0][deg_list[0].size()-1].emplace(gNameToID[deg_gene]);
                    }
//                    else if (regex_search(deg_gene, m, regex("(ENSMUSG)(.*)")))
                    else if (match_regex(& MM, &deg_gene[0]))
                    {
                        deg_list[1][deg_list[1].size()-1].emplace(gNameToID[deg_gene]);
                    }
                    else if (match_regex(& SB, &deg_gene[0]))
                    {
                        deg_list[2][deg_list[2].size()-1].emplace(gNameToID[deg_gene]);
                    }
                }
            }
            else
            {
                fprintf(stderr, "cannot open deg file\n");
            }
            DEGFILE.close();
            fprintf(stderr, "deg\t%s\tsize\t%d\n", line.c_str(), deg_count);
        }
    }
    else
    {
        fprintf(stderr, "cannot open deg list file\n");
    }
    DEGLIST.close();
}

int main(int argc, char * argv[])
{
    if (argc != 8)
    {
        printf("usage: degEnrich InputDir DEGList OrthGroups hbUniverse mmUniverse sbUniverse\n");
        exit(1);
    }
    /* read in universe gene sets */
    int arg_counter = 4;
    for (int i = 0; i<3;i++)
    {
        readUniverse(argv[arg_counter++]);
    }
    permute_round = atoi(argv[arg_counter++]);
    for (int i = 0; i<3;i++)
    {
        fprintf(stderr, "%d universe size\t%lu\n", i, universe_gene_set[i].size());
    }
    // read in Orth groups
    readOrthGroup(argv[3]);
    // read in DEG lists from all spes
    readDEGfile(argv[1], argv[2]);
    // generate permute universe and count orth edges for each pair of degs
    enrichInPermutUniverse();
    // find rank of orig pair of degs
    rankOfDGE();
    return 0;
}
