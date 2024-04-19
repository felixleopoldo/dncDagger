#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/tiernan_all_cycles.hpp>
#include <boost/graph/topological_sort.hpp>
#include <strings.h>
#include <string.h> // needed?
#include "auxiliary.h"

typedef boost::adjacency_matrix<boost::undirectedS> Graph;
typedef boost::adjacency_matrix<boost::directedS> DiGraph;


/* A struct for a subcomponent.
 */
struct SubComp
{
    vector<int> nodes;
    vector<int> suborder;    
    size_t max_n_particles = 0;
    size_t tot_n_particles = 0;
    vector<double> scores;
    vector<vector<bool>> opt_adjmat;
    vector<vector<bool>> subadjmat;
    double score = 0;
};

/* A struct for an isolated component.
 */
struct IsoComp
{
    vector<int> nodes;
    vector<int> subcompids;
    vector<double> scores;
    vector<SubComp> subcomps;
    bool connected = false;
    bool no_cycles = false;
    vector<vector<bool>> comp_dep;
    size_t max_n_particles = 0;
    size_t tot_n_particles = 0;
    double tot_order_to_dag_time = 0;
    double tot_order_time = 0;
};

struct IsoComps 
{
    vector<IsoComp> iso_comps;
    size_t max_n_particles = 0;
    size_t tot_n_particles = 0;
    double tot_order_to_dag_time = 0;
    double tot_order_time = 0;

};

vector<size_t> merged_neig_cycles(const vector<vector<bool>> & compdep);
vector<int> dnc(OrderScoring &scoring, vector<vector<bool>> &h_min, vector<vector<bool>> &h_max);
void printIsoComps(IsoComps &iso_comps);
bool restructure_components(IsoComp & iso_comp, OrderScoring & scoring);
void subcomponents_update(IsoComp &iso_comp, OrderScoring & scoring);
vector<int> concatenate_subcomponents(const IsoComp & iso_comp, OrderScoring & scoring);
vector<vector<bool>> subcomponents_dependence(const IsoComp & iso_comp, OrderScoring & scoring);
IsoComps structure_components(Graph &G_H_min, Graph &G_H_max);
pair<vector<vector<bool>>, vector<vector<bool>>> get_diff_matrices(vector<Rcpp::IntegerVector> rowmaps_backwards,vector<Rcpp::IntegerVector> rowmaps_forward, vector<vector<vector<double>>> scoretable, vector<vector<int>> potential_parents);
pair<vector<vector<bool>>, vector<vector<bool>>> get_diff_matrices(OrderScoring & scores);