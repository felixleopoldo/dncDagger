#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_traits.hpp>

typedef boost::adjacency_matrix<boost::directedS> Graph;


/* A struct for a subcomponent.
 */
struct SubComp
{
    vector<int> nodes;
    vector<int> suborder;
    vector<double> scores;
    vector<vector<int>> opt_adjmat;
    vector<vector<int>> subadjmat;
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
    int max_n_particles = 0;
    int tot_n_particles = 0;
    double tot_order_to_dag_time = 0;
    double tot_order_time = 0;
};

struct IsoComps 
{
    vector<IsoComp> iso_comps;
    int max_n_particles = 0;
    int tot_n_particles = 0;
    double tot_order_to_dag_time = 0;
    double tot_order_time = 0;
};


vector<vector<bool>> dnc(OrderScoring &scoring, vector<vector<bool>> &h_min, vector<vector<bool>> &h_max);

IsoComps structure_components(Graph &G_H_min, Graph &G_H_max);
