/* -------------------------------------------------------
AppPotts_eng class header for abnormal grain growth
--
Read-Shockley implementation developed by Efrain Hernandez-Rivera (2017--2018)
US Army Research Laboratory
--
THIS SOFTWARE IS MADE AVAILABLE ON AN "AS IS" BASIS
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, NEITHER
EXPRESSED OR IMPLIED
------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(potts/eng,AppPotts_eng)

#else

#ifndef SPK_APP_POTTS_ENG_H
#define SPK_APP_POTTS_ENG_H

#include <map>
#include <vector>
#include <unordered_map>
#include <string>
#include <set>
#include "app_potts.h"

namespace SPPARKS_NS {

class AppPotts_eng : public AppPotts {
  public:

    AppPotts_eng(class SPPARKS *, int, char **);
    ~AppPotts_eng();
    void input_app(char *, int, char **);
    void init_app();
    void grow_app();
    void site_event_rejection(int, class RandomPark *);

    void inclination_calculation(int i, int si, std::vector<double> & vector_lab_frame, int center);
    void compute_normal_vector_2DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center);
    int int_mod(int base, int divide);
    void vector_average(std::vector<double> vector_lab_frame_1, std::vector<double> vector_lab_frame_2, std::vector<double> & average_normal_lab_frame);
    void compute_normal_vector_3DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center);

    void convert_inclination_refercrystal(int i, std::vector<double> & average_normal_lab_frame, std::vector<double> & inclination);
    void vector_to_quternion(std::vector<double> & normals, std::vector<double> & reference_normals, std::vector<double> & q);

    void euler2quat(std::vector<double> eulerAngle, std::vector<double> & quternion_result);
    void mat2quat(const double O[3][3], double q[4]);
    void rotate_normals_byquat2mat(std::vector<double> q, std::vector<double> from, std::vector<double> & to);
    void symmat(double ***sym);
    void quat_mult(const double qi[4], const double qj[4], double q[4]);
    int in_cubic_fz(const double m_axis[3]);
    void quaternions_fz(const std::vector<double> qi, const std::vector<double> qj, std::vector<double> & q_result);

    // new neighbor function
    void get_nearest_neighbor(int i, int neighbor_size);
    void get_boarder_place(int center, int & i_dis, int & j_dis, int & k_dis);

    double compute_fully_energy_2d(std::vector<double>  misorientation, std::vector<double> inclination);
    double compute_fully_energy_5d(std::vector<double>  misorientation, std::vector<double> inclination);
    void energy_ave_split(int i, int nei, std::vector<double> & inclination);
    void misorientation_calculation(int i, int nei, std::vector<double> & misorientation);
    void read_linear_vector_matrix();
    double site_energy(int);





  protected:
    double *phi1,*phi2,*Phi; //pointer to 3 rotation angles
    int *spin;
    int dimension;  //simulation dimension
    double thetam; //High-low angle divider
    double Jij; //Interaction energy
    int Osym; //Symmetry Operator flag
    double **symquat; //Symmetry Operator in quaternion space
    //smooth algorithm parameters
    int interval, neighbor_length;
    std::string smthAlgo;
    int nx,ny,nz;
    // The calculated normal vector of site i or nei based on i's grain in lab frame
    std::vector<std::vector<double>> vector_i_lab_frame;
    // The calculated normal vector of site i or nei based on nei's grain in lab frame
    std::vector<std::vector<double>> vector_nei_lab_frame;
    // The averaged normal vector between site i and nei based on i's grain in lab frame
    std::vector<double> average_normal_i_lab_frame;
    // The averaged normal vector between site i and nei based on nei's grain in lab frame
    std::vector<double> average_normal_nei_lab_frame;

    std::vector<double> reference_axis = {0, 0, 0};

    // misorientation and inclination impact factor
    double m_impact_factor;
    double i_impact_factor;

    // inclination energy
    double E_delta;
    double E_m;
    double current_gamma = 0.0;

    // triple energy
    std::set<int> nei_set;
    std::string triple_energy;

    // vector matrix for L2
    std::vector<std::vector<double>> linear_vector_matrix_i;
    std::vector<std::vector<double>> linear_vector_matrix_j;
    std::vector<std::vector<std::vector<double>>> linear_vector_matrix_i_3d;
    std::vector<std::vector<std::vector<double>>> linear_vector_matrix_j_3d;
    std::vector<std::vector<std::vector<double>>> linear_vector_matrix_k_3d;

    // map to store misoreintation
    std::map<std::pair<int,int>, std::vector<double>> misorientation_map;

    // New neighbor functions
    int numneigh_near;
    std::vector<int> neighbor_near;


  private:
   double pfraction;
   int multi,nthresh;
};

}

#endif
#endif
