/* -------------------------------------------------------
AppPotts_eng class source for abnormal grain growth
--
Read-Shockley implementation developed by Efrain Hernandez-Rivera (2017--2018)
US Army Research Laboratory
--
THIS SOFTWARE IS MADE AVAILABLE ON AN "AS IS" BASIS
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, NEITHER
EXPRESSED OR IMPLIED
------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "domain.h"
#include "math.h"
#include "app.h"
#include "app_potts_eng.h"
#include "random_park.h"
#include "comm_lattice.h"
#include "error.h"
#include <fstream>
#include <iostream>
#include <type_traits>
#include <list>
using namespace SPPARKS_NS;

#define MY_PI 3.14159265358979323846 // pi
#define MY_2PI 6.28318530717958647692 // 2pi

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    out << "{";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
        if (i != last)
            out << ", ";
    }
    out << "}";
    return out;
}

/* ---------------------------------------------------- */

// I think this part is used for initialization
AppPotts_eng::AppPotts_eng(SPPARKS *spk, int narg, char **arg) :
AppPotts(spk,narg,arg)
{
  ninteger = 1;
  //1 double array per Euler angle
  ndouble = 3;
  // delpropensity = 6;
  // delevent = 0;
  // add the extra arrays
  recreate_arrays();
  // only error check for this class, not derived classes
  if (strcmp(arg[0],"potts/eng") == 0 && narg < 2)
  error->all(FLERR,"Illegal app_style command");
  //cutoff misorientation angle
  thetam=25.0/180.0*MY_PI;
  //interaction (interfacial) energy
  Jij=1.0;
  //Mobility parameters
  // nmob=4.0; bmob=5.0;
  nspins=atoi(arg[1]);
  //Symmetry operator
  Osym=24;
  if (narg == 3)
  Osym=atoi(arg[2]);
  //Inclination default parameters
  triple_energy = "ave";
  // smoothing algorithm iteration times
  interval = 3;
  neighbor_length = interval+1;
  reference_axis[0] = 1;
  // misorientation and inclination impact factor
  m_impact_factor = 1;
  i_impact_factor = 1;
}

/* -------------------------------------------------------
Destructor
------------------------------------------------------- */
AppPotts_eng::~AppPotts_eng()
{
  //free up memory from quaternion symmetry operator
  for (int i = 0; i<Osym; i++)
  delete[] symquat[i];
  delete[] symquat;
}

/* -------------------------------------------------------
Initialize before each run
check validity of site values
------------------------------------------------------- */
void AppPotts_eng::init_app()
{

  // the coordination of site i: xyz[i][0] - x; xyz[i][1] - y; xyz[i][2] - z;
  // xyz = app->xyz;

  int flag = 0;
  //Check angles are within corresponding range
  for (int i = 0; i < nlocal; i++) {
    if (phi1[i] < 0 || phi1[i] >= MY_2PI){
      fprintf(screen, "phi1 %d\n",i);
      flag = 1;
    }
    if (phi2[i] < 0 || phi2[i] >= MY_2PI) {
      fprintf(screen, "phi2 %d\n",i);
      flag = 1;
    }
    if (Phi[i] < 0 || Phi[i] >= MY_PI) {
      fprintf(screen, "Phi %d\n",i);
      flag = 1;
    }
    if (spin[i] < 1 || spin[i] > nspins+1) {
      fprintf(screen, "Spin %d\n",i);
      flag = 1;
    }
  }

  //Initialize symmetry operator as quaternion vectors
  //Osym = 24 (cubic), 12 (hexagonal)
  symmat(&symquat);

  comm->all();

  if (logfile)
    fprintf(logfile," Pairs misorientation map created\n");
  if (screen && me==0)
    fprintf(screen," Pairs misorientation map created\n");

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall)
    error->all(FLERR,"One or more sites have invalid values");

  // initial the parameters and matrix we need
  nx = ceil(domain->boxxhi);
  ny = ceil(domain->boxyhi);
  nz = ceil(domain->boxzhi);
  dimension = domain->dimension;

  if (dimension == 2) {
    linear_vector_matrix_i.assign(2*interval+3, std::vector<double>(2*interval+3, 0));
    linear_vector_matrix_j.assign(2*interval+3, std::vector<double>(2*interval+3, 0));

    // New neareast neighbor functions
    numneigh_near = 8;
    neighbor_near.assign(numneigh_near, 0);
  }
  else if (dimension == 3) {
    linear_vector_matrix_i_3d.assign(2*interval+3, std::vector<std::vector<double>>(2*interval+3, std::vector<double>(2*interval+3, 0)));
    linear_vector_matrix_j_3d.assign(2*interval+3, std::vector<std::vector<double>>(2*interval+3, std::vector<double>(2*interval+3, 0)));
    linear_vector_matrix_k_3d.assign(2*interval+3, std::vector<std::vector<double>>(2*interval+3, std::vector<double>(2*interval+3, 0)));

    // New neareast neighbor functions
    numneigh_near = 26;
    neighbor_near.assign(numneigh_near, 0);
  }

  delete [] sites; // include all object(like neighbors) the i site can flip
  delete [] unique;
  sites = new int[1 + numneigh_near];
  unique = new int[1 + numneigh_near];

  dt_sweep = 1.0/numneigh_near;

  // Inclination vector in lab and self frame assignment
  vector_i_lab_frame.assign(2, std::vector<double>(3,0));
  vector_nei_lab_frame.assign(2, std::vector<double>(3,0));
  average_normal_i_lab_frame.assign(3, 0.0);
  average_normal_nei_lab_frame.assign(3, 0.0);

  std::ifstream vectorMatrixFile("vectorMatrix/vectorMatrix_" +std::to_string(interval)+ ".txt");
  if (vectorMatrixFile.peek() == std::ifstream::traits_type::eof()) {
    // Save vectorMatrixFile list to file
    error->all(FLERR,"vectorMatrix files not exist.\n");
  }
  //Read in vectorMatrixFile
  else {
    read_linear_vector_matrix();
  }
  vectorMatrixFile.close();

  // create a output file to store the vector information
  std::ifstream vectorFile("Vector.txt");
  vectorFile.close();

}

void AppPotts_eng::read_linear_vector_matrix() {
  int ma_len = 2*interval+3;
  std::ifstream vectorMatrixFile("vectorMatrix/vectorMatrix_" +std::to_string(interval)+ ".txt");
  if (dimension==2) {
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        vectorMatrixFile >> linear_vector_matrix_i[i][j];
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        vectorMatrixFile >> linear_vector_matrix_j[i][j];
  }
  else if (dimension==3) {
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        for (int k=0; k<ma_len; k++)
          vectorMatrixFile >> linear_vector_matrix_i_3d[i][j][k];
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        for (int k=0; k<ma_len; k++)
          vectorMatrixFile >> linear_vector_matrix_j_3d[i][j][k];
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        for (int k=0; k<ma_len; k++)
          vectorMatrixFile >> linear_vector_matrix_k_3d[i][j][k];
  }
}

/* -------------------------------------------------------
Set site value ptrs each time iarray/darray are
reallocated
------------------------------------------------------- */
void AppPotts_eng::grow_app()
{
  // set pointers
  // to define these, use command
  // create_sites box iN and set iN
  spin = iarray[0];
  phi1 = darray[0];
  Phi = darray[1];
  phi2 = darray[2];
}

/* -------------------------------------------------------
User defined optional parameters
------------------------------------------------------- */
void AppPotts_eng::input_app(char *command, int narg, char **arg)
{
  //Redefine mobility parameters (n,b)
  if (strcmp(command,"mobility") == 0) {

  }
  //Cutoff angle for Read-Shockley
  else if (strcmp(command,"cutoff") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal cutoff angle command\n");
    thetam=fabs(atof(arg[0]))/180.0*MY_PI;
    if (thetam>MY_2PI)
      error->all(FLERR,"Cutoff angle must be defined in "
        "terms of degrees (0,360)\n");

    if (logfile)
      fprintf(logfile," Low-to-high angle cutoff reset "
        "to %s deg\n",arg[0]);
    if (screen && me==0)
      fprintf(screen," Low-to-high angle cutoff reset "
        "to %s deg\n",arg[0]);
  }
  //Potts interfacial energy scaler
  else if (strcmp(command,"energy_scaling") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal scaling energy command\n");
    Jij=atof(arg[0]);
    if (Jij<0)
      error->all(FLERR,"Illegal energy value (>0)\n");
    if (logfile)
      fprintf(logfile," PMC energy scaling by %g.\n",Jij);
    if (screen && me==0)
      fprintf(screen," PMC energy scaling by %g.\n",Jij);
  }
  else if (strcmp(command,"incParams") == 0) {
    if (narg != 4){
      error->all(FLERR,"Illegal incParams flag: requires "
        "two arguments, type double maximum inclination energy and type int number of discretized energy bins\n "
        "(e.g. incParams 1 5)\n");
    }
    else{
      E_delta = atof(arg[2]);
      E_m = atof(arg[3]);
    }
  }
  else if (strcmp(command,"interval") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal smoothing algorithm iteration command\n");
    else {
      interval=atof(arg[0]);
      neighbor_length = interval+1;
      smthAlgo = arg[1];
    }
    if (logfile)
      fprintf(logfile," We use smoothing algorithm %s with %d iteration.\n",smthAlgo.c_str(), interval);
    if (screen && me==0)
      fprintf(screen," We use smoothing algorithm %s with %d iteration.\n",smthAlgo.c_str(), interval);
  }
  else if (strcmp(command,"reference_axis") == 0) {
    if (narg<3)
      error->all(FLERR,"Illegal reference axis\n");
    reference_axis[0]=atof(arg[0]);
    reference_axis[1]=atof(arg[1]);
    reference_axis[2]=atof(arg[2]);
    if (logfile)
      fprintf(logfile," The reference axis is {%f, %f, %f}.\n", reference_axis[0], reference_axis[1], reference_axis[2]);
    if (screen && me==0)
      fprintf(screen," The reference axis is {%f, %f, %f}.\n", reference_axis[0], reference_axis[1], reference_axis[2]);
  }
  else if (strcmp(command,"TJ_energy_type") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal triple junction energy type\n");
    triple_energy=arg[0];
    if (logfile)
      fprintf(logfile," The constant triple energy is %s.\n", triple_energy.c_str());
    if (screen && me==0)
      fprintf(screen," The constant triple energy is %s.\n", triple_energy.c_str());
  }
  else if (strcmp(command,"impact_factor") == 0) {
    if (narg<2)
      error->all(FLERR,"Illegal misorientation and inclination impact factors\n");
    m_impact_factor=atof(arg[0]);
    i_impact_factor=atof(arg[1]);
    if (logfile)
      fprintf(logfile," The misorientation and inclination impact factors are %f and %f.\n", m_impact_factor, i_impact_factor);
    if (screen && me==0)
      fprintf(screen," The misorientation and inclination impact factors are %f and %f.\n", m_impact_factor, i_impact_factor);
  }
  else
    error->all(FLERR,"Input command not recognized by app\n");
}

void AppPotts_eng::get_nearest_neighbor(int i, int neighbor_size) {

  if (dimension == 2) {
    neighbor_near[0] = neighbor[i][2*neighbor_size*(neighbor_size+2)];
    neighbor_near[1] = neighbor[i][2*neighbor_size*(neighbor_size+2)+1];
    neighbor_near[2] = neighbor[i][2*neighbor_size*(neighbor_size+2)+2];
    neighbor_near[3] = neighbor[i][2*neighbor_size*neighbor_size+6*neighbor_size+3];
    neighbor_near[4] = neighbor[i][2*neighbor_size*neighbor_size+6*neighbor_size+4];
    neighbor_near[5] = neighbor[i][2*neighbor_size*neighbor_size+8*neighbor_size+5];
    neighbor_near[6] = neighbor[i][2*neighbor_size*neighbor_size+8*neighbor_size+6];
    neighbor_near[7] = neighbor[i][2*neighbor_size*neighbor_size+8*neighbor_size+7];
  }
  else if (dimension == 3) {
    neighbor_near[0] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+13*neighbor_size];
    neighbor_near[1] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+13*neighbor_size+1];
    neighbor_near[2] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+13*neighbor_size+2];
    neighbor_near[3] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+15*neighbor_size+3];
    neighbor_near[4] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+15*neighbor_size+4];
    neighbor_near[5] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+15*neighbor_size+5];
    neighbor_near[6] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+17*neighbor_size+6];
    neighbor_near[7] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+17*neighbor_size+7];
    neighbor_near[8] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+17*neighbor_size+8];

    neighbor_near[9] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+25*neighbor_size+9];
    neighbor_near[10] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+25*neighbor_size+10];
    neighbor_near[11] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+25*neighbor_size+11];
    neighbor_near[12] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+27*neighbor_size+12];
    neighbor_near[13] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+27*neighbor_size+13];
    neighbor_near[14] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+29*neighbor_size+14];
    neighbor_near[15] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+29*neighbor_size+15];
    neighbor_near[16] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+29*neighbor_size+16];

    neighbor_near[17] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+37*neighbor_size+17];
    neighbor_near[18] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+37*neighbor_size+18];
    neighbor_near[19] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+37*neighbor_size+19];
    neighbor_near[20] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+39*neighbor_size+20];
    neighbor_near[21] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+39*neighbor_size+21];
    neighbor_near[22] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+39*neighbor_size+22];
    neighbor_near[23] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+41*neighbor_size+23];
    neighbor_near[24] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+41*neighbor_size+24];
    neighbor_near[25] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+41*neighbor_size+25];
  }


}

/* -------------------------------------------------------
Compute Hamiltonian of site
------------------------------------------------------- */
double AppPotts_eng::site_energy(int i)
{
  if (spin[i] > nspins) return 0.0;
  int nei, nei_num = 0;
  double eng = 0.0, vec_len = 0.0;
  nei_set.clear();
  get_nearest_neighbor(i, interval);

  // If the site is not on the boundary
  for (int j = 0; j < numneigh_near; j++) {
    nei=neighbor_near[j];
    // Get neighboring grain ID
    if (spin[i] == spin[nei]) continue;
    nei_set.insert(spin[nei]);
    // Get num of neighboring site
    nei_num += 1;
  }
  if (!nei_set.empty()) {
    inclination_calculation(i, i, vector_i_lab_frame[0], i);
  }

  // energy
  // Boundary energy
  if (nei_set.size() == 1)
    for (int j = 0; j < numneigh_near; j++) {
      double tmp_eng=0;
      nei=neighbor_near[j];

      if (spin[i] == spin[nei]) continue;

      // Calculate misoreintation
      std::vector<double> misorientation = {0,0,0,0};
      misorientation_calculation(i, nei, misorientation);

      // add the inclination energy into total site energy
      std::vector<double> inclination = {0,0,0};
      energy_ave_split(i, nei, inclination);

      // Get energy
      tmp_eng += compute_fully_energy_5d(misorientation, inclination);
      eng+=tmp_eng;

    }
  else if (nei_set.size() > 1) {
    if (triple_energy == "min") {
      // Triple depend on lowest energy boundary site
      double min_eng = pow(3.0, dimension);
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // add the misorientation and inclination energy into total site energy
        std::vector<double> misorientation = {0,0,0,0};
        misorientation_calculation(i, nei, misorientation);
        std::vector<double> inclination = {0,0,0};
        energy_ave_split(i, nei, inclination);
        // Energy
        tmp_eng += compute_fully_energy_5d(misorientation, inclination);
        if (min_eng > tmp_eng) min_eng = tmp_eng;
      }
      eng += nei_num * min_eng;
    }
    else if (triple_energy == "max") {
      // Triple depend on highest energy boundary site
      double max_eng = 0;
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // add the misorientation and inclination energy into total site energy
        std::vector<double> misorientation = {0,0,0,0};
        misorientation_calculation(i, nei, misorientation);
        std::vector<double> inclination = {0,0,0};
        energy_ave_split(i, nei, inclination);
        // Energy
        tmp_eng += compute_fully_energy_5d(misorientation, inclination);
        if (max_eng < tmp_eng) max_eng = tmp_eng;
      }
      eng += nei_num * max_eng;
    }
    else if (triple_energy == "ave") {
      // New AVE TJ energy type
      double sum_eng = 0;
      std::vector<std::vector<double>> other_grain_inclination(nei_set.size(), std::vector<double>(3, 0)); // store inclination for other grains except for spin[i]
      std::vector<int> other_grain_num(nei_set.size(), 0); // store num of sites for other grains except for spin[i]
      for (int j = 0; j < numneigh_near; j++) {
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // Get the inclination for nei
        std::vector<double> inclination = {0,0,0};
        energy_ave_split(i, nei, inclination);
        // Store the orientation and inclination and num of sites into corresponding grain sequence
        std::set<int>::iterator iter;
        int k = 0;
        for (iter = nei_set.begin(); iter != nei_set.end(); ++iter) {
          if (spin[nei] == *iter) {
            for (int m=0; m < 3; m++) {
              other_grain_inclination[k][m] += vector_nei_lab_frame[1][m];
            }
            other_grain_num[k] += 1;
          }
          k++;
        }
      }
      // Get true orientation and average inclination for each grain
      for (int k = 0; k < other_grain_inclination.size(); k++) {
        for (int m=0; m < other_grain_inclination[k].size(); m++) {
          other_grain_inclination[k][m] = other_grain_inclination[k][m] / other_grain_num[k];
        }
        vec_len = sqrt(other_grain_inclination[k][0]*other_grain_inclination[k][0]+other_grain_inclination[k][1]*other_grain_inclination[k][1]+other_grain_inclination[k][2]*other_grain_inclination[k][2]);
        if (vec_len > 1e-5) {
          double ff1 = 1/vec_len;
          other_grain_inclination[k][0] = other_grain_inclination[k][0]*ff1;
          other_grain_inclination[k][1] = other_grain_inclination[k][1]*ff1;
          other_grain_inclination[k][2] = other_grain_inclination[k][2]*ff1;
        }
        else {
          other_grain_inclination[k] = {0,0,1};
        }
      }
      // Get Ave energy
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng=0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // Calculate misoreintation
        std::vector<double> misorientation = {0,0,0,0};
        misorientation_calculation(i, nei, misorientation);
        // Couple with all average inclination
        for (int k = 0; k < other_grain_inclination.size(); k++) {
          tmp_eng += compute_fully_energy_5d(misorientation, other_grain_inclination[k]);
        }
        tmp_eng += compute_fully_energy_5d(misorientation, vector_i_lab_frame[0]); // Take the i inclination into account
        sum_eng += tmp_eng / (nei_set.size() + 1);
      }
      eng += sum_eng;
    }
    else if (triple_energy == "consMin") {
      // Triple depend on constant boundary site (minimal energy)
      double global_min_eng = 0;
      global_min_eng = 1 - E_delta;
      eng += nei_num * global_min_eng;
    }
    else if (triple_energy == "consMax") {
      // Triple depend on constant boundary site (minimal energy)
      double global_max_eng = 0;
      global_max_eng = 1 + E_delta;
      eng += nei_num * global_max_eng;
    }
    else if (triple_energy == "consTest") {
      // Triple depend on constant boundary site (minimal energy)
      double global_constest_eng = 0.093;
      eng += nei_num * global_constest_eng;
    }
    else if (triple_energy == "old") {
      // Triple depend on constant boundary site (minimal energy)

      for (int j = 0; j < numneigh_near; j++) {
        double old_TJ_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // add the misorientation and inclination energy into total site energy
        std::vector<double> misorientation = {0,0,0,0};
        misorientation_calculation(i, nei, misorientation);
        //  inclination
        std::vector<double> inclination = {0,0,0};
        energy_ave_split(i, nei, inclination);
        // Energy
        old_TJ_eng = compute_fully_energy_5d(misorientation, inclination);
        eng += old_TJ_eng;
      }
    }
    else {
      error->all(FLERR,"Please use the correct TJ energy approach\n");
    }
  }
  current_gamma = 1.0;

  return Jij*eng;
}

void AppPotts_eng::energy_ave_split(int i, int nei, std::vector<double> & inclination) {
  // Calculatie inclination

  // Calculate the normal vector of i's grain in lab frame
  inclination_calculation(i, nei, vector_i_lab_frame[1], i);
  vector_average(vector_i_lab_frame[0], vector_i_lab_frame[1], average_normal_i_lab_frame); // Get the average normals for sites in two sides with lab frame
  // Calculate the normal vector of nei's grain in lab frame
  inclination_calculation(nei, nei, vector_nei_lab_frame[1], i);
  inclination_calculation(nei, i, vector_nei_lab_frame[0], i);
  vector_average(vector_nei_lab_frame[1], vector_nei_lab_frame[0], average_normal_nei_lab_frame);
  vector_average(average_normal_i_lab_frame, average_normal_nei_lab_frame, inclination); // final inclination
}

void AppPotts_eng::misorientation_calculation(int i, int nei, std::vector<double> & misorientation) {
  // Calculate misorientation from i and nei

  // misorientation pair with min-max
  int smin = MIN(spin[i],spin[nei]);
  int smax = MAX(spin[i],spin[nei]);
  std::pair <int, int> spins = std::make_pair(smin,smax);
  std::map<std::pair<int,int>, std::vector<double>>::iterator map_iter;
  map_iter = misorientation_map.find(spins);
  // Add misorientation dict if exist i-nei
  if (map_iter != misorientation_map.end()) {
    misorientation = map_iter->second;
  }
  else { // Add misorientation dict if not exist
    // vector initialization
    std::vector<double> orientation_i_ea = {phi1[i], Phi[i], phi2[i]};
    std::vector<double> orientation_nei_ea = {phi1[nei], Phi[nei], phi2[nei]};
    std::vector<double> orientation_i_qua = {0,0,0,0};
    std::vector<double> orientation_nei_qua = {0,0,0,0};
    // calculation misorientation
    euler2quat(orientation_i_ea, orientation_i_qua);
    euler2quat(orientation_nei_ea, orientation_nei_qua);
    if (spin[i]<spin[nei]) quaternions_fz(orientation_i_qua, orientation_nei_qua, misorientation);
    else quaternions_fz(orientation_nei_qua, orientation_i_qua, misorientation);
    // add into dict
    misorientation_map[spins] = misorientation;
  }

}

void AppPotts_eng::euler2quat(std::vector<double> eulerAngle, std::vector<double> & quternion_result)
{
    //Convert grain Euler angles to quaternion vector
    quternion_result[0]=cos(eulerAngle[1]/2.)*cos((eulerAngle[0]+eulerAngle[2])/2.);
    quternion_result[1]=sin(eulerAngle[1]/2.)*cos((eulerAngle[0]-eulerAngle[2])/2.);
    quternion_result[2]=sin(eulerAngle[1]/2.)*sin((eulerAngle[0]-eulerAngle[2])/2.);
    quternion_result[3]=cos(eulerAngle[1]/2.)*sin((eulerAngle[0]+eulerAngle[2])/2.);

}

/* -------------------------------------------------------
Convert symmetry matrix to quaternion form
------------------------------------------------------- */
void AppPotts_eng::mat2quat(const double O[3][3], double q[4])
{
    double q4 = 0;
    if ( (1 + O[0][0] + O[1][1] + O[2][2]) > 0) {
        q4 = sqrt(1 + O[0][0] + O[1][1] + O[2][2])/2;
        q[0] = q4;
        q[1] = (O[2][1] - O[1][2])/(4*q4);
        q[2] = (O[0][2] - O[2][0])/(4*q4);
        q[3] = (O[1][0] - O[0][1])/(4*q4);
    }
    else if ( (1 + O[0][0] - O[1][1] - O[2][2]) > 0) {
        q4 = sqrt(1 + O[0][0] - O[1][1] - O[2][2])/2;
        q[0] = (O[2][1] - O[1][2])/(4*q4);
        q[1] = q4;
        q[2] = (O[1][0] + O[0][1])/(4*q4);
        q[3] = (O[0][2] + O[2][0])/(4*q4);
    }
    else if ( (1 - O[0][0] + O[1][1] - O[2][2]) > 0) {
        q4 = sqrt(1 - O[0][0] + O[1][1] - O[2][2])/2;
        q[0] = (O[0][2] - O[2][0])/(4*q4);
        q[1] = (O[1][0] + O[0][1])/(4*q4);
        q[2] = q4;
        q[3] = (O[2][1] + O[1][2])/(4*q4);
    }
    else if ( (1 - O[0][0] - O[1][1] + O[2][2]) > 0) {
        q4 = sqrt(1 - O[0][0] - O[1][1] + O[2][2])/2;
        q[0] = (O[1][0] - O[0][1])/(4*q4);
        q[1] = (O[0][2] + O[2][0])/(4*q4);
        q[2] = (O[2][1] + O[1][2])/(4*q4);
        q[3] = q4;
    }
}
void AppPotts_eng::rotate_normals_byquat2mat(std::vector<double> q, std::vector<double> from, std::vector<double> & to)
{
    std::vector<std::vector<double>> O = {{0,0,0},{0,0,0},{0,0,0}};

    O[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    O[0][1] = 2*(q[1]*q[2] - q[0]*q[3]);
    O[0][2] = 2*(q[0]*q[2] + q[1]*q[3]);
    O[1][0] = 2*(q[1]*q[2] + q[0]*q[3]);
    O[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    O[1][2] = 2*(-q[0]*q[1] + q[2]*q[3]);
    O[2][0] = 2*(-q[0]*q[2] + q[1]*q[3]);
    O[2][1] = 2*(q[0]*q[1] + q[2]*q[3]);
    O[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];

    to[0] = O[0][0] * from[0] + O[0][1] * from[1] + O[0][2] * from[2];
    to[1] = O[1][0] * from[0] + O[1][1] * from[1] + O[1][2] * from[2];
    to[2] = O[2][0] * from[0] + O[2][1] * from[1] + O[2][2] * from[2];

}

/* -------------------------------------------------------
Define the symmetry operator
------------------------------------------------------- */
void AppPotts_eng::symmat(double ***sym)
{
    //grow by number of symmetric operators
    (*sym) = new double*[Osym];

    //grow for symmetry quaternion vectors
    for (int o=0; o<Osym; o++)
        (*sym)[o] = new double[4];

    //buffer for quaternion
    double q[4];

    if (Osym == 24) {
        //cubic symmetry
        double SYM[24][3][3] =
                { {{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}},
                  {{ 1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}},
                  {{ 1, 0, 0}, { 0, 0,-1}, { 0, 1, 0}},
                  {{ 1, 0, 0}, { 0, 0, 1}, { 0,-1, 0}},
                  {{-1, 0, 0}, { 0, 1, 0}, { 0, 0,-1}},
                  {{-1, 0, 0}, { 0,-1, 0}, { 0, 0, 1}},
                  {{-1, 0, 0}, { 0, 0,-1}, { 0,-1, 0}},
                  {{-1, 0, 0}, { 0, 0, 1}, { 0, 1, 0}},
                  {{ 0, 1, 0}, {-1, 0, 0}, { 0, 0, 1}},
                  {{ 0, 1, 0}, { 0, 0,-1}, {-1, 0, 0}},
                  {{ 0, 1, 0}, { 1, 0, 0}, { 0, 0,-1}},
                  {{ 0, 1, 0}, { 0, 0, 1}, { 1, 0, 0}},
                  {{ 0,-1, 0}, { 1, 0, 0}, { 0, 0, 1}},
                  {{ 0,-1, 0}, { 0, 0,-1}, { 1, 0, 0}},
                  {{ 0,-1, 0}, {-1, 0, 0}, { 0, 0,-1}},
                  {{ 0,-1, 0}, { 0, 0, 1}, {-1, 0, 0}},
                  {{ 0, 0, 1}, { 0, 1, 0}, {-1, 0, 0}},
                  {{ 0, 0, 1}, { 1, 0, 0}, { 0, 1, 0}},
                  {{ 0, 0, 1}, { 0,-1, 0}, { 1, 0, 0}},
                  {{ 0, 0, 1}, {-1, 0, 0}, { 0,-1, 0}},
                  {{ 0, 0,-1}, { 0, 1, 0}, { 1, 0, 0}},
                  {{ 0, 0,-1}, {-1, 0, 0}, { 0, 1, 0}},
                  {{ 0, 0,-1}, { 0,-1, 0}, {-1, 0, 0}},
                  {{ 0, 0,-1}, { 1, 0, 0}, { 0,-1, 0}} };

        //initialize global operator
        for (int o=0; o<Osym; o++) {
            mat2quat(SYM[o],q);
            for (int i=0; i<4; i++)
                (*sym)[o][i]=q[i];
        }
    }
    else if (Osym == 12) {
        //hexagonal symmetry
        double a = sqrt(3)/2;
        double SYM[12][3][3] =
                { {{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}},
                  {{-0.5, a, 0}, { -a,-0.5, 0}, { 0, 0, 1}},
                  {{-0.5, -a, 0}, { a,-0.5, 0}, { 0, 0, 1}},
                  {{ 0.5, a, 0}, { -a, 0.5, 0}, { 0, 0, 1}},
                  {{ -1, 0, 0}, { 0, -1, 0}, { 0, 0, 1}},
                  {{ 0.5, -a, 0}, { a, 0.5, 0}, { 0, 0, 1}},
                  {{-0.5, -a, 0}, { -a, 0.5, 0}, { 0, 0, -1}},
                  {{ 1, 0, 0}, { 0, -1, 0}, { 0, 0, -1}},
                  {{-0.5, a, 0}, { a, 0.5, 0}, { 0, 0, -1}},
                  {{ 0.5, a, 0}, { a,-0.5, 0}, { 0, 0, -1}},
                  {{ -1, 0, 0}, { 0, 1, 0}, { 0, 0, -1}},
                  {{ 0.5, -a, 0}, { -a,-0.5, 0}, { 0, 0, -1}} };

        //initialize global operator
        for (int o=0; o<Osym; o++) {
            mat2quat(SYM[o],q);
            for (int i=0; i<4; i++)
                (*sym)[o][i]=q[i];
        }
    }
}

void AppPotts_eng::quat_mult(const double qi[4], const double qj[4], double q[4])
{
  //Hamilton multiplication/product
  //multiplying quaternions and update
  q[0] = qi[0]*qj[0] - qi[1]*qj[1] - qi[2]*qj[2] - qi[3]*qj[3];
  q[1] = qi[0]*qj[1] + qi[1]*qj[0] + qi[2]*qj[3] - qi[3]*qj[2];
  q[2] = qi[0]*qj[2] - qi[1]*qj[3] + qi[2]*qj[0] + qi[3]*qj[1];
  q[3] = qi[0]*qj[3] + qi[1]*qj[2] - qi[2]*qj[1] + qi[3]*qj[0];
}

int AppPotts_eng::in_cubic_fz(const double m_axis[3])
{
  // check if the misoreintation axis in fundamental zone

  // three core axis
  double axis_a[3] = {1,0,0}, axis_b[3] = {1/sqrt(2),1/sqrt(2),0}, axis_c[3] = {1/sqrt(3),1/sqrt(3),1/sqrt(3)};
  // if in fundamental zone
  bool judgement_0 = ((axis_a[1]*axis_b[2] - axis_a[2]*axis_b[1])*m_axis[0] + (axis_a[2]*axis_b[0] - axis_a[0]*axis_b[2])*m_axis[1] + (axis_a[0]*axis_b[1] - axis_a[1]*axis_b[0])*m_axis[2] >= 0);
  bool judgement_1 = ((axis_b[1]*axis_c[2] - axis_b[2]*axis_c[1])*m_axis[0] + (axis_b[2]*axis_c[0] - axis_b[0]*axis_c[2])*m_axis[1] + (axis_b[0]*axis_c[1] - axis_b[1]*axis_c[0])*m_axis[2] >= 0);
  bool judgement_2 = ((axis_c[1]*axis_a[2] - axis_c[2]*axis_a[1])*m_axis[0] + (axis_c[2]*axis_a[0] - axis_c[0]*axis_a[2])*m_axis[1] + (axis_c[0]*axis_a[1] - axis_c[1]*axis_a[0])*m_axis[2] >= 0);

  return int(judgement_0 * judgement_1 * judgement_2);
}

void AppPotts_eng::quaternions_fz(const std::vector<double> qi, const std::vector<double> qj, std::vector<double> & q_result)
{
  double miso0, misom=MY_2PI, m_axis_base;
  q_result.assign(4, 0);

  double q[4], qib[4], qjb[4], qmin[4]={0,0,0,0}, m_axis[3]={1,0,0};
  double tmp_q[4]={0,0,0,0};
  double qi_array[4], qj_array[4];
  qi_array[0]=qi[0]; qi_array[1]=qi[1]; qi_array[2]=qi[2]; qi_array[3]=qi[3];
  qj_array[0]=qj[0]; qj_array[1]=qj[1]; qj_array[2]=qj[2]; qj_array[3]=qj[3];
  for (int o1=0; o1<Osym; o1++) {
    for (int o2=0; o2<Osym; o2++) {
      quat_mult(symquat[o1],qi_array,qib);
      quat_mult(symquat[o2],qj_array,qjb);

      //j-grain conjugate quaternion
      qjb[1]=-qjb[1];
      qjb[2]=-qjb[2];
      qjb[3]=-qjb[3];
      quat_mult(qib,qjb,q);
      // misorientation axis in fundamental zone
      m_axis_base = sqrt(1-q[0]*q[0]);
      if (m_axis_base>0) {
        m_axis[0] = q[1]/m_axis_base;
        m_axis[1] = q[2]/m_axis_base;
        m_axis[2] = q[3]/m_axis_base;
      }
      // check if the q and inverse of q in fundamental zone
      if (!in_cubic_fz(m_axis)) {
        q[1]=-q[1];
        q[2]=-q[2];
        q[3]=-q[3];
        m_axis[0] = -m_axis[0];
        m_axis[1] = -m_axis[1];
        m_axis[2] = -m_axis[2];
        if (!in_cubic_fz(m_axis)) continue;
      }
      // if in fundamental zone, find the q with minimal m-angle
      miso0 = 2*acos(round(q[0]*1e5)/1e5);
      if (miso0 > MY_PI) miso0 = miso0-MY_2PI;
      if (fabs(miso0) < misom) {
        misom=fabs(miso0);
        qmin[0]=q[0]; qmin[1]=q[1]; qmin[2]=q[2]; qmin[3]=q[3];
      }

    }
  }

  q_result[0] = qmin[0];
  q_result[1] = qmin[1];
  q_result[2] = qmin[2];
  q_result[3] = qmin[3];

  return;
}

/* -------------------------------------------------------
Calculate the rotation of quaternion version
between two vectors (from normals to reference_normals)
------------------------------------------------------- */
void AppPotts_eng::vector_to_quternion(std::vector<double> & normals, std::vector<double> & reference_normals, std::vector<double> & q)
{
  // Initialization
  q.assign(4, 0);
  double cos_theta;

  // normals cross product reference_normals to get the rotation axis
  q[1] = normals[1]*reference_normals[2] - normals[2]*reference_normals[1];
  q[2] = normals[2]*reference_normals[0] - normals[0]*reference_normals[2];
  q[3] = normals[0]*reference_normals[1] - normals[1]*reference_normals[0];



  // normals with reference_normals to get the rotation angle
  cos_theta = round((normals[0]*reference_normals[0] + normals[1]*reference_normals[1] + normals[2]*reference_normals[2])*1e5)/1e5;
  q[0] = sqrt((cos_theta+1)/2);

  // get the quaternion by the rotation angle and axis
  q[1] = q[1]*sqrt((1-cos_theta)/2);
  q[2] = q[2]*sqrt((1-cos_theta)/2);
  q[3] = q[3]*sqrt((1-cos_theta)/2);

}


/* -------------------------------------------------------
Calculate the middle point vector by
averaging the two site vectors
------------------------------------------------------- */
void AppPotts_eng::vector_average(std::vector<double> vector_lab_frame_1, std::vector<double> vector_lab_frame_2, std::vector<double> & average_normal_lab_frame)
{
  average_normal_lab_frame.assign(3,0);
  average_normal_lab_frame[0] = vector_lab_frame_1[0] - vector_lab_frame_2[0];
  average_normal_lab_frame[1] = vector_lab_frame_1[1] - vector_lab_frame_2[1];
  average_normal_lab_frame[2] = vector_lab_frame_1[2] - vector_lab_frame_2[2];

  double normal_len = sqrt(average_normal_lab_frame[0]*average_normal_lab_frame[0]+average_normal_lab_frame[1]*average_normal_lab_frame[1]+average_normal_lab_frame[2]*average_normal_lab_frame[2]);

  if (normal_len < 1e-6) {
    average_normal_lab_frame[0] = vector_lab_frame_1[0];
    average_normal_lab_frame[1] = vector_lab_frame_1[1];
    average_normal_lab_frame[2] = vector_lab_frame_1[2];
    return;
  }
  average_normal_lab_frame[0] /= normal_len;
  average_normal_lab_frame[1] /= normal_len;
  average_normal_lab_frame[2] /= normal_len;

  return;
}

/* -------------------------------------------------------
Calculate the inclination on si by the crystal frame of i
------------------------------------------------------- */
void AppPotts_eng::inclination_calculation(int i, int si, std::vector<double> & vector_lab_frame, int center)
{
  // output the inclination vector by matrix method
  vector_lab_frame.assign(3,0);

  switch (dimension) {
    case 2:
      compute_normal_vector_2DLinear_matrix(i, si, vector_lab_frame, center);
      break;
    case 3:
      compute_normal_vector_3DLinear_matrix(i, si, vector_lab_frame, center);
      break;
  }
}

/* -------------------------------------------------------
Calculate the inclination from the vector
------------------------------------------------------- */
void AppPotts_eng::convert_inclination_refercrystal(int i, std::vector<double> & average_normal_lab_frame, std::vector<double> & inclination)
{
  inclination.assign(3,0);
  inclination[0] = 1.0 * (cos(phi1[i])*cos(phi2[i])-sin(phi1[i])*sin(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[0] +\
                              (-cos(phi1[i])*sin(phi2[i])-sin(phi1[i])*cos(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[1] +\
                              (sin(phi1[i])*sin(Phi[i])) * average_normal_lab_frame[2];
  inclination[1] = 1.0 * (sin(phi1[i])*cos(phi2[i])+cos(phi1[i])*sin(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[0] +\
                              (-sin(phi1[i])*sin(phi2[i])+cos(phi1[i])*cos(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[1] +\
                              (-cos(phi1[i])*sin(Phi[i])) * average_normal_lab_frame[2];
  inclination[2] = 1.0 * (sin(phi2[i])*sin(Phi[i])) * average_normal_lab_frame[0] +\
                              (cos(phi2[i])*sin(Phi[i])) * average_normal_lab_frame[1] +\
                              (cos(Phi[i])) * average_normal_lab_frame[2];
  return;
}

/* -------------------------------------------------------
Calculate the inclination energy at a voxel
------------------------------------------------------- */
double AppPotts_eng::compute_fully_energy_2d(std::vector<double>  misorientation, std::vector<double> inclination)
{

  double inc_energy;
  double max_theta = 180.0/180.0*MY_PI;
  double m_axis_base = sqrt(1-misorientation[0]*misorientation[0]);
  std::vector<double> misorientation_axis = {1,1,1};
  if (m_axis_base != 0) {
    misorientation_axis = {misorientation[1]/m_axis_base,\
                           misorientation[2]/m_axis_base,\
                           misorientation[3]/m_axis_base};
  }
  double misorientation_angle = 2*acos(abs(misorientation[0]));
  double cos_inclination = round((misorientation_axis[0]*inclination[0] + misorientation_axis[1]*inclination[1] + misorientation_axis[2]*inclination[2])*1e5)/1e5;
  double inclination_angle = acos(abs(cos_inclination));

  // twist boundary (two axis are same)
  double energy_twist = 0.2 + 0.8 * misorientation_angle / max_theta * (1 - log(misorientation_angle / max_theta));
  // tile boundary (two axis are perpendeicular)
  double energy_tilt = 0.2 + 0.3 * misorientation_angle / max_theta * (1 - log(misorientation_angle / max_theta));

  // Get energy by inclination angle
  if (inclination_angle > MY_PI/2) inclination_angle = MY_PI - inclination_angle;
  inc_energy = energy_twist + inclination_angle / (MY_PI/2) * (energy_tilt - energy_twist);

  return inc_energy;
}
double AppPotts_eng::compute_fully_energy_5d(std::vector<double>  misorientation, std::vector<double> inclination)
{

  double inc_energy;
  double max_theta = 65.0/180.0*MY_PI;// 10.0/180.0*MY_PI;
  double m_axis_base = sqrt(1-misorientation[0]*misorientation[0]);
  std::vector<double> misorientation_axis = {1,0,0}; // if m_axis_base == 0
  if (m_axis_base != 0) {
    misorientation_axis = {misorientation[1]/m_axis_base,
                           misorientation[2]/m_axis_base,
                           misorientation[3]/m_axis_base};
  }
  double misorientation_angle = 2*acos(abs(misorientation[0]));  // [0, pi]
  double cos_inclination = round((misorientation_axis[0]*inclination[0] + misorientation_axis[1]*inclination[1] + misorientation_axis[2]*inclination[2])*1e5)/1e5;
  double inclination_angle = acos(abs(cos_inclination));  // [0, pi]
  if (misorientation_angle > max_theta) misorientation_angle = max_theta;

  // The two angular parameters for misorientation axis (m-axis - either orientation)
  double m_polar_angle = acos(round(misorientation_axis[2]*1e5)/1e5); // [0, pi]
  double m_azimuth_angle = atan2(misorientation_axis[1], misorientation_axis[0]) + MY_PI; // [0, 2*pi]

  // The two angular parameters for inclination axis (i-axis - m-axis)
  std::vector<double> rotation_quaternion_from_maxis_to_zaxis = {0,0,0,0};
  std::vector<double> inclination_axis_refer_misorientation_axis = {0,0,0};
  std::vector<double> z_axis = {0,0,1};
  vector_to_quternion(misorientation_axis, z_axis, rotation_quaternion_from_maxis_to_zaxis);
  rotate_normals_byquat2mat(rotation_quaternion_from_maxis_to_zaxis, inclination, inclination_axis_refer_misorientation_axis);
  double i_polar_angle = acos(round(inclination_axis_refer_misorientation_axis[2]*1e5)/1e5); // [0, pi]
  double i_azimuth_angle = atan2(inclination_axis_refer_misorientation_axis[1], inclination_axis_refer_misorientation_axis[0]) + MY_PI; // [0, 2*pi]

  // misorientation axis effect
  // double m_axis_effect = abs(cos(m_azimuth_angle)) * m_polar_angle / MY_PI;
  double m_axis_effect = pow(abs(cos(m_azimuth_angle/2)),0.4) + pow(abs(cos(m_polar_angle)),0.4); // new m_axis_effect with circle impact
  // double m_axis_effect = (abs(m_azimuth_angle-MY_PI) < MY_PI/18 && abs(m_polar_angle - MY_PI/2) < MY_PI/36)? 0.0: 1.0; // new m_axis_effect with circle impact
  if (m_axis_effect > 1) m_axis_effect = 1;
  // twist boundary (two axis are same) and tile boundary (two axis are perpendeicular)
  double energy_twist = 0.7 * (1.0-m_impact_factor); // if misorientation == 0
  double energy_tilt = 0.3 * (1.0-m_impact_factor); // if misorientation == 0
  if (misorientation_angle != 0) {
    energy_twist = 0.7 * (m_axis_effect * misorientation_angle / max_theta * (1 - log(misorientation_angle / max_theta))*m_impact_factor + (1.0-m_impact_factor));
    energy_tilt = 0.3 * (m_axis_effect *  misorientation_angle / max_theta * (1 - log(misorientation_angle / max_theta))*m_impact_factor + (1.0-m_impact_factor));
  }

  // Get energy by inclination angle
  if (i_polar_angle > MY_PI/2) i_polar_angle = MY_PI - i_polar_angle;
  inc_energy = 0.3 + abs(cos(i_azimuth_angle)) * (energy_twist + i_polar_angle / (MY_PI/2) * (energy_tilt - energy_twist)) * i_impact_factor + energy_twist * (1.0-i_impact_factor);
  // double inc_axis_effect = ((abs(i_azimuth_angle-MY_PI*0.5) < MY_PI/18) || (abs(i_azimuth_angle-MY_PI*1.5) < MY_PI/18)) && (abs(i_polar_angle-MY_PI/4)< MY_PI/36)? 0.0: 1.0;
  // double inc_axis_effect = ((abs(i_azimuth_angle-MY_PI) < MY_PI/18) || (abs(i_azimuth_angle) < MY_PI/18)) && (abs(i_polar_angle-MY_PI/4)< MY_PI/36)? 0.0: 1.0;
  // inc_energy = 0.3 + inc_axis_effect * energy_twist * i_impact_factor + energy_twist * (1.0-i_impact_factor);

  return inc_energy;
}

void AppPotts_eng::get_boarder_place(int center, int & i_dis, int & j_dis, int & k_dis) {
  double x_min, y_min, z_min;
  double x_ind, y_ind, z_ind;
  x_min = 1.0*(domain->boxxhi - domain->boxxlo) / domain->procgrid[0];
  y_min = 1.0*(domain->boxyhi - domain->boxylo) / domain->procgrid[1];
  z_min = 1.0*(domain->boxzhi - domain->boxzlo) / domain->procgrid[2];
  i_dis = 0;
  j_dis = 0;
  k_dis = 0;

  if (dimension == 2) {
    x_ind = xyz[center][0];
    y_ind = xyz[center][1];
    while (x_ind>=x_min) x_ind = x_ind - x_min;
    while (y_ind>=y_min) y_ind = y_ind - y_min;

    if (x_ind < 1) j_dis = -1;
    else if (x_min - x_ind <= 1) j_dis = 1;
    if (y_ind < 1) i_dis = -1;
    else if (y_min - y_ind <= 1) i_dis = 1;

    if (i_dis==0 && j_dis==0) std::cout << "*** SPPARKS cannot find the correct boarder place 1 *** " << xyz[center][0] << " " << xyz[center][1] << std::endl;

  }
  else if (dimension == 3) {
    x_ind = xyz[center][0];
    y_ind = xyz[center][1];
    z_ind = xyz[center][2];
    while (x_ind>=x_min) x_ind = x_ind - x_min;
    while (y_ind>=y_min) y_ind = y_ind - y_min;
    while (z_ind>=z_min) z_ind = z_ind - z_min;

    if (x_ind < 1) j_dis = -1;
    else if (x_min - x_ind <= 1) j_dis = 1;
    if (y_ind < 1) i_dis = -1;
    else if (y_min - y_ind <= 1) i_dis = 1;
    if (z_ind < 1) k_dis = -1;
    else if (z_min - z_ind <= 1) k_dis = 1;

    if (i_dis==0 && j_dis==0 && k_dis==0) std::cout << "*** SPPARKS cannot find the correct boarder place *** " << xyz[center][0] << " " << xyz[center][1] << " " << xyz[center][2] << std::endl;

  }

  return;
}

void AppPotts_eng::compute_normal_vector_2DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center)
{
  // Get the normals at si, all spin == spin[i] is 1

  int sval, local_x, local_y;
  int j_dis = 0, i_dis = 0, k_dis = 0;
  double inc_x=0, inc_y=0, vec_len;

  if (numneigh[si] == (2*interval+3)*(2*interval+3)-1){
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j;
      if (ni < ((2*interval+3)*(2*interval+3)-1)/2) {
        norm_i = ni / (2*interval+3);
        norm_j = ni - (2*interval+3)*norm_i;
      }
      else {
        norm_i = (ni+1) / (2*interval+3);
        norm_j = (ni+1) - (2*interval+3)*norm_i;
      }

      inc_x += factor_i*linear_vector_matrix_i[norm_i][norm_j];
      inc_y += factor_i*linear_vector_matrix_j[norm_i][norm_j];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  } else if (numneigh[si] == (2*interval+3)*(2*interval+2)-1) {
    get_boarder_place(center, i_dis, j_dis, k_dis);
    if (i_dis*j_dis!=0) {
      if (xyz[si][1]==xyz[center][1]) i_dis = 0;
      if (xyz[si][0]==xyz[center][0]) j_dis = 0;
    }
    if (i_dis*j_dis!=0) {
      i_dis = 0; j_dis = 0; k_dis = 0;
      get_boarder_place(si, i_dis, j_dis, k_dis);
      i_dis = -i_dis; j_dis = -j_dis; k_dis = -k_dis;
    }
    if (i_dis*j_dis!=0) std::cout <<"Unexpected bool issue " << "ij_dis: " << i_dis<<" "<<j_dis << " si: " << xyz[si][0] << " " << xyz[si][1] << " center: " << xyz[center][0] << " " << xyz[center][1] << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j;
      if (ni < (interval+1-(i_dis<0?1:0))*(2*interval+3-abs(j_dis))+interval+1-(j_dis<0?1:0)) {
        norm_i = ni / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = ni - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_i = (ni+1) / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = (ni+1) - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }

      inc_x += factor_i*linear_vector_matrix_i[norm_i][norm_j];
      inc_y += factor_i*linear_vector_matrix_j[norm_i][norm_j];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  } else if (numneigh[si] == (2*interval+2)*(2*interval+2)-1) {
    get_boarder_place(si, i_dis, j_dis, k_dis);
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j;
      if (ni < (interval+1-(i_dis<0?1:0))*(2*interval+3-abs(j_dis))+interval+1-(j_dis<0?1:0)) {
        norm_i = ni / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = ni - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_i = (ni+1) / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = (ni+1) - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }

      inc_x += factor_i*linear_vector_matrix_i[norm_i][norm_j];
      inc_y += factor_i*linear_vector_matrix_j[norm_i][norm_j];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  }

  vec_len = sqrt(inc_x*inc_x+inc_y*inc_y);
  if (vec_len > 1e-5) {
    vector_lab_frame[0] = -inc_y/vec_len;
    vector_lab_frame[1] = -inc_x/vec_len;
    vector_lab_frame[2] = 0;
  }

  return;
}

void AppPotts_eng::compute_normal_vector_3DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center)
{
  // Get the normals at si, all spin == spin[i] is 1

  int sval, local_x, local_y, local_z;;
  int j_dis = 0, i_dis = 0, k_dis = 0;
  int j_dis_si = 0, i_dis_si = 0, k_dis_si = 0;
  double inc_x=0, inc_y=0, inc_z=0, vec_len;

  if (numneigh[si] == (2*interval+3)*(2*interval+3)*(2*interval+3)-1){
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < ((2*interval+3)*(2*interval+3)*(2*interval+3)-1)/2) {
        norm_k = ni / ((2*interval+3)*(2*interval+3));
        norm_i = (ni - ((2*interval+3)*(2*interval+3))*norm_k) / (2*interval+3);
        norm_j = ni - ((2*interval+3)*(2*interval+3))*norm_k - (2*interval+3) * norm_i;
      }
      else {
        norm_k = (ni+1) / ((2*interval+3)*(2*interval+3));
        norm_i = ((ni+1) - ((2*interval+3)*(2*interval+3))*norm_k) / (2*interval+3);
        norm_j = (ni+1) - ((2*interval+3)*(2*interval+3))*norm_k - (2*interval+3) * norm_i;
      }
      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_k_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  } else if (numneigh[si] == (2*interval+3)*(2*interval+3)*(2*interval+2)-1) {
    get_boarder_place(center, i_dis, j_dis, k_dis);
    if (!(bool(i_dis)^bool(j_dis)^bool(k_dis)) || bool(i_dis*j_dis*k_dis)) {
      if (xyz[si][1]==xyz[center][1]) i_dis = 0;
      if (xyz[si][0]==xyz[center][0]) j_dis = 0;
      if (xyz[si][2]==xyz[center][2]) k_dis = 0;
    }
    if (!(bool(i_dis)^bool(j_dis)^bool(k_dis)) || bool(i_dis*j_dis*k_dis)) {
      i_dis_si = 0; j_dis_si = 0; k_dis_si = 0;
      get_boarder_place(si, i_dis_si, j_dis_si, k_dis_si);
      i_dis_si = -i_dis_si; j_dis_si = -j_dis_si; k_dis_si = -k_dis_si;
      i_dis = i_dis * (i_dis == i_dis_si);
      j_dis = j_dis * (j_dis == j_dis_si);
      k_dis = k_dis * (k_dis == k_dis_si);
    }
    if (!(bool(i_dis)^bool(j_dis)^bool(k_dis)) || bool(i_dis*j_dis*k_dis)) std::cout << "Unexpected bool issue!!" << "jik_dis: " << j_dis<<" "<<i_dis <<" "<<k_dis<<" jik_dis_si: "<<j_dis_si<<" "<<i_dis_si<<" "<<k_dis_si<< " si: " << xyz[si][0] << " " << xyz[si][1] <<" "<<xyz[si][2] << " center: " << xyz[center][0] << " " << xyz[center][1]<<" "<<xyz[center][2] << " expected: "<<numneigh[si] << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(interval+1-(k_dis<0?1:0))+
      (2*interval+3-abs(j_dis))*(interval+1-(i_dis<0?1:0))+interval+1-(j_dis<0?1:0)) {
        norm_k = ni / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = (ni -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = ni - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_k = (ni+1) / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = ((ni+1) -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = (ni+1) - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_k_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  } else if (numneigh[si] == (2*interval+3)*(2*interval+2)*(2*interval+2)-1) {
    get_boarder_place(center, i_dis, j_dis, k_dis);
    if (!((bool(i_dis) && bool(j_dis)) ^ (bool(j_dis) && bool(k_dis)) ^ (bool(i_dis) && bool(k_dis))) || bool(i_dis*j_dis*k_dis)) {
      if (xyz[si][1]==xyz[center][1]) i_dis = 0;
      if (xyz[si][0]==xyz[center][0]) j_dis = 0;
      if (xyz[si][2]==xyz[center][2]) k_dis = 0;
    }
    if (!((bool(i_dis) && bool(j_dis)) ^ (bool(j_dis) && bool(k_dis)) ^ (bool(i_dis) && bool(k_dis))) || bool(i_dis*j_dis*k_dis)) {
      i_dis_si = 0; j_dis_si = 0; k_dis_si = 0;
      get_boarder_place(si, i_dis_si, j_dis_si, k_dis_si);
      i_dis_si = -i_dis_si; j_dis_si = -j_dis_si; k_dis_si = -k_dis_si;
      i_dis = i_dis * (i_dis == i_dis_si);
      j_dis = j_dis * (j_dis == j_dis_si);
      k_dis = k_dis * (k_dis == k_dis_si);
    }
    if (!((bool(i_dis) && bool(j_dis)) ^ (bool(j_dis) && bool(k_dis)) ^ (bool(i_dis) && bool(k_dis))) || bool(i_dis*j_dis*k_dis)) std::cout << "Unexpected bool issue!!" << "jik_dis: " << j_dis<<" "<<i_dis <<" "<<k_dis<<" jik_dis_si: "<<j_dis_si<<" "<<i_dis_si<<" "<<k_dis_si<< " si: " << xyz[si][0] << " " << xyz[si][1] <<" "<<xyz[si][2] << " center: " << xyz[center][0] << " " << xyz[center][1]<<" "<<xyz[center][2] << " expected: "<<numneigh[si] << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(interval+1-(k_dis<0?1:0))+
      (2*interval+3-abs(j_dis))*(interval+1-(i_dis<0?1:0))+interval+1-(j_dis<0?1:0)) {
        norm_k = ni / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = (ni -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = ni - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_k = (ni+1) / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = ((ni+1) -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = (ni+1) - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }

      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_k_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  } else if (numneigh[si] == (2*interval+2)*(2*interval+2)*(2*interval+2)-1) {
    get_boarder_place(center, i_dis, j_dis, k_dis);
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(interval+1-(k_dis<0?1:0))+
      (2*interval+3-abs(j_dis))*(interval+1-(i_dis<0?1:0))+interval+1-(j_dis<0?1:0)) {
        norm_k = ni / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = (ni -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = ni - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_k = (ni+1) / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = ((ni+1) -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = (ni+1) - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_k_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  }


  vec_len = sqrt(inc_x*inc_x+inc_y*inc_y+inc_z*inc_z);
  if (vec_len > 1e-5) {
    double ff1 = 1/vec_len;
    vector_lab_frame[0] = -inc_y*ff1;
    vector_lab_frame[1] = -inc_x*ff1;
    vector_lab_frame[2] = -inc_z*ff1;
  }

  return;
}

int AppPotts_eng::int_mod(int base, int divide) {
  return int(base - floor(double(base)/divide)*divide);
}

/* -------------------------------------------------------
rKMC method
perform a site event with no null bin rejection
flip to random neighbor spin without null bin
------------------------------------------------------- */
void AppPotts_eng::site_event_rejection(int i, RandomPark *random)
{
  // no events for a pinned site
  if (spin[i] > nspins) return;

  int oldstate=spin[i];
  double iphi[3]={phi1[i],Phi[i],phi2[i]};

  // events = spin flips to neighboring site different than self

  int j,nei;
  int nevent = 0;
  float tmp_error =0.0;
  get_nearest_neighbor(i, interval);

  //Nearest-neighbor sampling
  for (j = 0; j < numneigh_near; j++) {
    nei=neighbor_near[j];
    if (spin[i]==spin[nei])
      continue;
    sites[nevent++]=nei;
  }

  if (nevent == 0) return;
  int iran = (int) (nevent*random->uniform());
  if (iran >= nevent) iran = nevent-1;

  double einitial = site_energy(i);

  spin[i] = spin[sites[iran]];
  phi1[i] = phi1[sites[iran]];
  phi2[i] = phi2[sites[iran]];
  Phi[i] = Phi[sites[iran]];

  double efinal = site_energy(i);

  // now the p0 is 1 in aniso
  double p0=1;
  // double p0=current_gamma;

  // accept or reject via Boltzmann criterion
  if (efinal <= einitial) {
    if (random->uniform() < p0) {
    }
    else {
      spin[i] = oldstate;
      phi1[i] = iphi[0];
      phi2[i] = iphi[2];
      Phi[i] = iphi[1];
    }
  }
  else if (temperature == 0.0) {
    spin[i] = oldstate;
    phi1[i] = iphi[0];
    phi2[i] = iphi[2];
    Phi[i] = iphi[1];
  }
  else if (random->uniform() > p0*exp((einitial-efinal)*t_inverse/current_gamma)) {
    spin[i] = oldstate;
    phi1[i] = iphi[0];
    phi2[i] = iphi[2];
    Phi[i] = iphi[1];
  }
  else {
  }

  if (spin[i] != oldstate) {
    naccept++;
  }
}
