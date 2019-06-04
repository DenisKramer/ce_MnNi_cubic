#include <cstddef>
#include "casm/clex/Clexulator.hh"



/****** CLEXULATOR CLASS FOR PRIM ******
{
  "basis" : [
    {
      "coordinate" : [ 0.000000000000, 0.000000000000, 0.000000000000 ],
      "occupant_dof" : [ "La" ]
    },
    {
      "coordinate" : [ 0.500000000000, 0.500000000000, 0.500000000000 ],
      "occupant_dof" : [ "Mn", "Ni" ]
    },
    {
      "coordinate" : [ 0.000000000000, 0.500000000000, 0.500000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.500000000000, 0.000000000000, 0.500000000000 ],
      "occupant_dof" : [ "O" ]
    },
    {
      "coordinate" : [ 0.500000000000, 0.500000000000, 0.000000000000 ],
      "occupant_dof" : [ "O" ]
    }
  ],
  "coordinate_mode" : "Fractional",
  "lattice_vectors" : [
    [ 3.940144706398, 0.000000000000, 0.000000000000 ],
    [ 0.000000000000, 3.940144706398, 0.000000000000 ],
    [ 0.000000000000, 0.000000000000, 3.940144706398 ]
  ],
  "title" : "La_MnNi_O3_cubic"
}**/


/// \brief Returns a Clexulator_impl::Base* owning a La_MnNi_O3_cubic_Clexulator
extern "C" CASM::Clexulator_impl::Base* make_La_MnNi_O3_cubic_Clexulator();

namespace CASM {

  class La_MnNi_O3_cubic_Clexulator : public Clexulator_impl::Base {

  public:

    La_MnNi_O3_cubic_Clexulator();

    ~La_MnNi_O3_cubic_Clexulator();

    /// \brief Clone the La_MnNi_O3_cubic_Clexulator
    std::unique_ptr<La_MnNi_O3_cubic_Clexulator> clone() const { 
      return std::unique_ptr<La_MnNi_O3_cubic_Clexulator>(_clone()); 
    }

    /// \brief Calculate contribution to global correlations from one unit cell
    void calc_global_corr_contribution(double *corr_begin) const override;

    /// \brief Calculate contribution to select global correlations from one unit cell
    void calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;

    /// \brief Calculate point correlations about basis site 'b_index'
    void calc_point_corr(int b_index, double *corr_begin) const override;

    /// \brief Calculate select point correlations about basis site 'b_index'
    void calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;

    /// \brief Calculate the change in point correlations due to changing an occupant
    void calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const override;

    /// \brief Calculate the change in select point correlations due to changing an occupant
    void calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const override;


  private:

    /// \brief Clone the Clexulator
    virtual La_MnNi_O3_cubic_Clexulator* _clone() const override {
      return new La_MnNi_O3_cubic_Clexulator(*this);
    }

    // typedef for method pointers
    typedef double (La_MnNi_O3_cubic_Clexulator::*BasisFuncPtr)() const;

    // typedef for method pointers
    typedef double (La_MnNi_O3_cubic_Clexulator::*DeltaBasisFuncPtr)(int, int) const;

    // array of pointers to member functions for calculating basis functions
    BasisFuncPtr m_orbit_func_list[35];

    // array of pointers to member functions for calculating flower functions
    BasisFuncPtr m_flower_func_lists[5][35];

    // array of pointers to member functions for calculating DELTA flower functions
    DeltaBasisFuncPtr m_delta_func_lists[5][35];

    // Occupation Function tables for basis sites in asymmetric unit 1:
    //   - basis site 1:
    double m_occ_func_1_0[2];

    // Occupation Function accessors for basis site 1:
    const double &occ_func_1_0(const int &nlist_ind)const{return m_occ_func_1_0[*(m_occ_ptr+*(m_nlist_ptr+nlist_ind))];}

    //default functions for basis function evaluation 
    double zero_func() const{ return 0.0;};
    double zero_func(int,int) const{ return 0.0;};

    double eval_bfunc_0_0_0() const;

    double eval_bfunc_1_0_0() const;

    double site_eval_at_1_bfunc_1_0_0() const;

    double delta_site_eval_at_1_bfunc_1_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_0_0() const;

    double site_eval_at_1_bfunc_2_0_0() const;

    double delta_site_eval_at_1_bfunc_2_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_1_0() const;

    double site_eval_at_1_bfunc_2_1_0() const;

    double delta_site_eval_at_1_bfunc_2_1_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_2_0() const;

    double site_eval_at_1_bfunc_2_2_0() const;

    double delta_site_eval_at_1_bfunc_2_2_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_3_0() const;

    double site_eval_at_1_bfunc_2_3_0() const;

    double delta_site_eval_at_1_bfunc_2_3_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_4_0() const;

    double site_eval_at_1_bfunc_2_4_0() const;

    double delta_site_eval_at_1_bfunc_2_4_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_5_0() const;

    double site_eval_at_1_bfunc_2_5_0() const;

    double delta_site_eval_at_1_bfunc_2_5_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_6_0() const;

    double site_eval_at_1_bfunc_2_6_0() const;

    double delta_site_eval_at_1_bfunc_2_6_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_7_0() const;

    double site_eval_at_1_bfunc_2_7_0() const;

    double delta_site_eval_at_1_bfunc_2_7_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_8_0() const;

    double site_eval_at_1_bfunc_2_8_0() const;

    double delta_site_eval_at_1_bfunc_2_8_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_9_0() const;

    double site_eval_at_1_bfunc_2_9_0() const;

    double delta_site_eval_at_1_bfunc_2_9_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_10_0() const;

    double site_eval_at_1_bfunc_2_10_0() const;

    double delta_site_eval_at_1_bfunc_2_10_0(int occ_i, int occ_f) const;

    double eval_bfunc_2_11_0() const;

    double site_eval_at_1_bfunc_2_11_0() const;

    double delta_site_eval_at_1_bfunc_2_11_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_0_0() const;

    double site_eval_at_1_bfunc_3_0_0() const;

    double delta_site_eval_at_1_bfunc_3_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_1_0() const;

    double site_eval_at_1_bfunc_3_1_0() const;

    double delta_site_eval_at_1_bfunc_3_1_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_2_0() const;

    double site_eval_at_1_bfunc_3_2_0() const;

    double delta_site_eval_at_1_bfunc_3_2_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_3_0() const;

    double site_eval_at_1_bfunc_3_3_0() const;

    double delta_site_eval_at_1_bfunc_3_3_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_4_0() const;

    double site_eval_at_1_bfunc_3_4_0() const;

    double delta_site_eval_at_1_bfunc_3_4_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_5_0() const;

    double site_eval_at_1_bfunc_3_5_0() const;

    double delta_site_eval_at_1_bfunc_3_5_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_6_0() const;

    double site_eval_at_1_bfunc_3_6_0() const;

    double delta_site_eval_at_1_bfunc_3_6_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_7_0() const;

    double site_eval_at_1_bfunc_3_7_0() const;

    double delta_site_eval_at_1_bfunc_3_7_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_8_0() const;

    double site_eval_at_1_bfunc_3_8_0() const;

    double delta_site_eval_at_1_bfunc_3_8_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_9_0() const;

    double site_eval_at_1_bfunc_3_9_0() const;

    double delta_site_eval_at_1_bfunc_3_9_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_10_0() const;

    double site_eval_at_1_bfunc_3_10_0() const;

    double delta_site_eval_at_1_bfunc_3_10_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_11_0() const;

    double site_eval_at_1_bfunc_3_11_0() const;

    double delta_site_eval_at_1_bfunc_3_11_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_12_0() const;

    double site_eval_at_1_bfunc_3_12_0() const;

    double delta_site_eval_at_1_bfunc_3_12_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_13_0() const;

    double site_eval_at_1_bfunc_3_13_0() const;

    double delta_site_eval_at_1_bfunc_3_13_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_14_0() const;

    double site_eval_at_1_bfunc_3_14_0() const;

    double delta_site_eval_at_1_bfunc_3_14_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_15_0() const;

    double site_eval_at_1_bfunc_3_15_0() const;

    double delta_site_eval_at_1_bfunc_3_15_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_16_0() const;

    double site_eval_at_1_bfunc_3_16_0() const;

    double delta_site_eval_at_1_bfunc_3_16_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_17_0() const;

    double site_eval_at_1_bfunc_3_17_0() const;

    double delta_site_eval_at_1_bfunc_3_17_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_18_0() const;

    double site_eval_at_1_bfunc_3_18_0() const;

    double delta_site_eval_at_1_bfunc_3_18_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_19_0() const;

    double site_eval_at_1_bfunc_3_19_0() const;

    double delta_site_eval_at_1_bfunc_3_19_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_20_0() const;

    double site_eval_at_1_bfunc_3_20_0() const;

    double delta_site_eval_at_1_bfunc_3_20_0(int occ_i, int occ_f) const;


  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  La_MnNi_O3_cubic_Clexulator::La_MnNi_O3_cubic_Clexulator() :
    Clexulator_impl::Base(179, 35) {
    m_occ_func_1_0[0] = -0.0000000000, m_occ_func_1_0[1] = 1.0000000000;

    m_orbit_func_list[0] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_0_0_0;
    m_orbit_func_list[1] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_1_0_0;
    m_orbit_func_list[2] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_0_0;
    m_orbit_func_list[3] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_1_0;
    m_orbit_func_list[4] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_2_0;
    m_orbit_func_list[5] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_3_0;
    m_orbit_func_list[6] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_4_0;
    m_orbit_func_list[7] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_5_0;
    m_orbit_func_list[8] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_6_0;
    m_orbit_func_list[9] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_7_0;
    m_orbit_func_list[10] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_8_0;
    m_orbit_func_list[11] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_9_0;
    m_orbit_func_list[12] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_10_0;
    m_orbit_func_list[13] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_11_0;
    m_orbit_func_list[14] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_0_0;
    m_orbit_func_list[15] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_1_0;
    m_orbit_func_list[16] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_2_0;
    m_orbit_func_list[17] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_3_0;
    m_orbit_func_list[18] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_4_0;
    m_orbit_func_list[19] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_5_0;
    m_orbit_func_list[20] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_6_0;
    m_orbit_func_list[21] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_7_0;
    m_orbit_func_list[22] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_8_0;
    m_orbit_func_list[23] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_9_0;
    m_orbit_func_list[24] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_10_0;
    m_orbit_func_list[25] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_11_0;
    m_orbit_func_list[26] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_12_0;
    m_orbit_func_list[27] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_13_0;
    m_orbit_func_list[28] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_14_0;
    m_orbit_func_list[29] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_15_0;
    m_orbit_func_list[30] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_16_0;
    m_orbit_func_list[31] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_17_0;
    m_orbit_func_list[32] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_18_0;
    m_orbit_func_list[33] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_19_0;
    m_orbit_func_list[34] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_20_0;


    m_flower_func_lists[0][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[0][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_flower_func_lists[1][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[1][1] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_1_0_0;
    m_flower_func_lists[1][2] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_0_0;
    m_flower_func_lists[1][3] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_1_0;
    m_flower_func_lists[1][4] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_2_0;
    m_flower_func_lists[1][5] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_3_0;
    m_flower_func_lists[1][6] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_4_0;
    m_flower_func_lists[1][7] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_5_0;
    m_flower_func_lists[1][8] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_6_0;
    m_flower_func_lists[1][9] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_7_0;
    m_flower_func_lists[1][10] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_8_0;
    m_flower_func_lists[1][11] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_9_0;
    m_flower_func_lists[1][12] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_10_0;
    m_flower_func_lists[1][13] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_11_0;
    m_flower_func_lists[1][14] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_0_0;
    m_flower_func_lists[1][15] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_1_0;
    m_flower_func_lists[1][16] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_2_0;
    m_flower_func_lists[1][17] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_3_0;
    m_flower_func_lists[1][18] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_4_0;
    m_flower_func_lists[1][19] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_5_0;
    m_flower_func_lists[1][20] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_6_0;
    m_flower_func_lists[1][21] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_7_0;
    m_flower_func_lists[1][22] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_8_0;
    m_flower_func_lists[1][23] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_9_0;
    m_flower_func_lists[1][24] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_10_0;
    m_flower_func_lists[1][25] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_11_0;
    m_flower_func_lists[1][26] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_12_0;
    m_flower_func_lists[1][27] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_13_0;
    m_flower_func_lists[1][28] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_14_0;
    m_flower_func_lists[1][29] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_15_0;
    m_flower_func_lists[1][30] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_16_0;
    m_flower_func_lists[1][31] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_17_0;
    m_flower_func_lists[1][32] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_18_0;
    m_flower_func_lists[1][33] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_19_0;
    m_flower_func_lists[1][34] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_20_0;


    m_flower_func_lists[2][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[2][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_flower_func_lists[3][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[3][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_flower_func_lists[4][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[4][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_delta_func_lists[0][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[0][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_delta_func_lists[1][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[1][1] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_1_0_0;
    m_delta_func_lists[1][2] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_0_0;
    m_delta_func_lists[1][3] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_1_0;
    m_delta_func_lists[1][4] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_2_0;
    m_delta_func_lists[1][5] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_3_0;
    m_delta_func_lists[1][6] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_4_0;
    m_delta_func_lists[1][7] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_5_0;
    m_delta_func_lists[1][8] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_6_0;
    m_delta_func_lists[1][9] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_7_0;
    m_delta_func_lists[1][10] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_8_0;
    m_delta_func_lists[1][11] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_9_0;
    m_delta_func_lists[1][12] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_10_0;
    m_delta_func_lists[1][13] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_11_0;
    m_delta_func_lists[1][14] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_0_0;
    m_delta_func_lists[1][15] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_1_0;
    m_delta_func_lists[1][16] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_2_0;
    m_delta_func_lists[1][17] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_3_0;
    m_delta_func_lists[1][18] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_4_0;
    m_delta_func_lists[1][19] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_5_0;
    m_delta_func_lists[1][20] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_6_0;
    m_delta_func_lists[1][21] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_7_0;
    m_delta_func_lists[1][22] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_8_0;
    m_delta_func_lists[1][23] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_9_0;
    m_delta_func_lists[1][24] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_10_0;
    m_delta_func_lists[1][25] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_11_0;
    m_delta_func_lists[1][26] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_12_0;
    m_delta_func_lists[1][27] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_13_0;
    m_delta_func_lists[1][28] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_14_0;
    m_delta_func_lists[1][29] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_15_0;
    m_delta_func_lists[1][30] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_16_0;
    m_delta_func_lists[1][31] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_17_0;
    m_delta_func_lists[1][32] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_18_0;
    m_delta_func_lists[1][33] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_19_0;
    m_delta_func_lists[1][34] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_20_0;


    m_delta_func_lists[2][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[2][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_delta_func_lists[3][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[3][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_delta_func_lists[4][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][1] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][2] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][3] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][4] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][5] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][6] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][7] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][8] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][9] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][10] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][11] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][12] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][13] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][14] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][15] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][16] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][17] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][18] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][19] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][20] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][21] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][22] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][23] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][24] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][25] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][26] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][27] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][28] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][29] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][30] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][31] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][32] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][33] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[4][34] = &La_MnNi_O3_cubic_Clexulator::zero_func;


    m_weight_matrix.row(0) << 1, 0, 0;
    m_weight_matrix.row(1) << 0, 1, 0;
    m_weight_matrix.row(2) << 0, 0, 1;

    m_neighborhood = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -3, -1, -1)},
      {UnitCellCoord(1, -3, -1, 0)},
      {UnitCellCoord(1, -3, -1, 1)},
      {UnitCellCoord(1, -3, 0, -1)},
      {UnitCellCoord(1, -3, 0, 0)},
      {UnitCellCoord(1, -3, 0, 1)},
      {UnitCellCoord(1, -3, 1, -1)},
      {UnitCellCoord(1, -3, 1, 0)},
      {UnitCellCoord(1, -3, 1, 1)},
      {UnitCellCoord(1, -2, -2, -2)},
      {UnitCellCoord(1, -2, -2, -1)},
      {UnitCellCoord(1, -2, -2, 0)},
      {UnitCellCoord(1, -2, -2, 1)},
      {UnitCellCoord(1, -2, -2, 2)},
      {UnitCellCoord(1, -2, -1, -2)},
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, -1, 2)},
      {UnitCellCoord(1, -2, 0, -2)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 0, 2)},
      {UnitCellCoord(1, -2, 1, -2)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -2, 1, 2)},
      {UnitCellCoord(1, -2, 2, -2)},
      {UnitCellCoord(1, -2, 2, -1)},
      {UnitCellCoord(1, -2, 2, 0)},
      {UnitCellCoord(1, -2, 2, 1)},
      {UnitCellCoord(1, -2, 2, 2)},
      {UnitCellCoord(1, -1, -3, -1)},
      {UnitCellCoord(1, -1, -3, 0)},
      {UnitCellCoord(1, -1, -3, 1)},
      {UnitCellCoord(1, -1, -2, -2)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -2, 2)},
      {UnitCellCoord(1, -1, -1, -3)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, -1, 3)},
      {UnitCellCoord(1, -1, 0, -3)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 0, 3)},
      {UnitCellCoord(1, -1, 1, -3)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 1, 3)},
      {UnitCellCoord(1, -1, 2, -2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, -1, 2, 2)},
      {UnitCellCoord(1, -1, 3, -1)},
      {UnitCellCoord(1, -1, 3, 0)},
      {UnitCellCoord(1, -1, 3, 1)},
      {UnitCellCoord(1, 0, -3, -1)},
      {UnitCellCoord(1, 0, -3, 0)},
      {UnitCellCoord(1, 0, -3, 1)},
      {UnitCellCoord(1, 0, -2, -2)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -2, 2)},
      {UnitCellCoord(1, 0, -1, -3)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, -1, 3)},
      {UnitCellCoord(1, 0, 0, -3)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 0, 3)},
      {UnitCellCoord(1, 0, 1, -3)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 1, 3)},
      {UnitCellCoord(1, 0, 2, -2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 0, 2, 2)},
      {UnitCellCoord(1, 0, 3, -1)},
      {UnitCellCoord(1, 0, 3, 0)},
      {UnitCellCoord(1, 0, 3, 1)},
      {UnitCellCoord(1, 1, -3, -1)},
      {UnitCellCoord(1, 1, -3, 0)},
      {UnitCellCoord(1, 1, -3, 1)},
      {UnitCellCoord(1, 1, -2, -2)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -2, 2)},
      {UnitCellCoord(1, 1, -1, -3)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, -1, 3)},
      {UnitCellCoord(1, 1, 0, -3)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 0, 3)},
      {UnitCellCoord(1, 1, 1, -3)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 1, 1)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 1, 3)},
      {UnitCellCoord(1, 1, 2, -2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 1, 2, 2)},
      {UnitCellCoord(1, 1, 3, -1)},
      {UnitCellCoord(1, 1, 3, 0)},
      {UnitCellCoord(1, 1, 3, 1)},
      {UnitCellCoord(1, 2, -2, -2)},
      {UnitCellCoord(1, 2, -2, -1)},
      {UnitCellCoord(1, 2, -2, 0)},
      {UnitCellCoord(1, 2, -2, 1)},
      {UnitCellCoord(1, 2, -2, 2)},
      {UnitCellCoord(1, 2, -1, -2)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, -1, 2)},
      {UnitCellCoord(1, 2, 0, -2)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 0)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 0, 2)},
      {UnitCellCoord(1, 2, 1, -2)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 0)},
      {UnitCellCoord(1, 2, 1, 1)},
      {UnitCellCoord(1, 2, 1, 2)},
      {UnitCellCoord(1, 2, 2, -2)},
      {UnitCellCoord(1, 2, 2, -1)},
      {UnitCellCoord(1, 2, 2, 0)},
      {UnitCellCoord(1, 2, 2, 1)},
      {UnitCellCoord(1, 2, 2, 2)},
      {UnitCellCoord(1, 3, -1, -1)},
      {UnitCellCoord(1, 3, -1, 0)},
      {UnitCellCoord(1, 3, -1, 1)},
      {UnitCellCoord(1, 3, 0, -1)},
      {UnitCellCoord(1, 3, 0, 0)},
      {UnitCellCoord(1, 3, 0, 1)},
      {UnitCellCoord(1, 3, 1, -1)},
      {UnitCellCoord(1, 3, 1, 0)},
      {UnitCellCoord(1, 3, 1, 1)}
    };


    m_orbit_neighborhood.resize(corr_size());
    m_orbit_neighborhood[0] = std::set<UnitCellCoord> {
    };

    m_orbit_neighborhood[1] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, 0, 0, 0)}
    };

    m_orbit_neighborhood[2] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 1, 0, 0)}
    };

    m_orbit_neighborhood[3] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, 0)}
    };

    m_orbit_neighborhood[4] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 1)}
    };

    m_orbit_neighborhood[5] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 2, 0, 0)}
    };

    m_orbit_neighborhood[6] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, 0)}
    };

    m_orbit_neighborhood[7] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[8] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -2, 0)},
      {UnitCellCoord(1, -2, 0, -2)},
      {UnitCellCoord(1, -2, 0, 2)},
      {UnitCellCoord(1, -2, 2, 0)},
      {UnitCellCoord(1, 0, -2, -2)},
      {UnitCellCoord(1, 0, -2, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 2, -2)},
      {UnitCellCoord(1, 0, 2, 2)},
      {UnitCellCoord(1, 2, -2, 0)},
      {UnitCellCoord(1, 2, 0, -2)},
      {UnitCellCoord(1, 2, 0, 2)},
      {UnitCellCoord(1, 2, 2, 0)}
    };

    m_orbit_neighborhood[9] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -2, -1)},
      {UnitCellCoord(1, -2, -2, 1)},
      {UnitCellCoord(1, -2, -1, -2)},
      {UnitCellCoord(1, -2, -1, 2)},
      {UnitCellCoord(1, -2, 1, -2)},
      {UnitCellCoord(1, -2, 1, 2)},
      {UnitCellCoord(1, -2, 2, -1)},
      {UnitCellCoord(1, -2, 2, 1)},
      {UnitCellCoord(1, -1, -2, -2)},
      {UnitCellCoord(1, -1, -2, 2)},
      {UnitCellCoord(1, -1, 2, -2)},
      {UnitCellCoord(1, -1, 2, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 1, -2, -2)},
      {UnitCellCoord(1, 1, -2, 2)},
      {UnitCellCoord(1, 1, 2, -2)},
      {UnitCellCoord(1, 1, 2, 2)},
      {UnitCellCoord(1, 2, -2, -1)},
      {UnitCellCoord(1, 2, -2, 1)},
      {UnitCellCoord(1, 2, -1, -2)},
      {UnitCellCoord(1, 2, -1, 2)},
      {UnitCellCoord(1, 2, 1, -2)},
      {UnitCellCoord(1, 2, 1, 2)},
      {UnitCellCoord(1, 2, 2, -1)},
      {UnitCellCoord(1, 2, 2, 1)}
    };

    m_orbit_neighborhood[10] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -3, 0, 0)},
      {UnitCellCoord(1, 0, -3, 0)},
      {UnitCellCoord(1, 0, 0, -3)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 3)},
      {UnitCellCoord(1, 0, 3, 0)},
      {UnitCellCoord(1, 3, 0, 0)}
    };

    m_orbit_neighborhood[11] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -3, -1, 0)},
      {UnitCellCoord(1, -3, 0, -1)},
      {UnitCellCoord(1, -3, 0, 1)},
      {UnitCellCoord(1, -3, 1, 0)},
      {UnitCellCoord(1, -1, -3, 0)},
      {UnitCellCoord(1, -1, 0, -3)},
      {UnitCellCoord(1, -1, 0, 3)},
      {UnitCellCoord(1, -1, 3, 0)},
      {UnitCellCoord(1, 0, -3, -1)},
      {UnitCellCoord(1, 0, -3, 1)},
      {UnitCellCoord(1, 0, -1, -3)},
      {UnitCellCoord(1, 0, -1, 3)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -3)},
      {UnitCellCoord(1, 0, 1, 3)},
      {UnitCellCoord(1, 0, 3, -1)},
      {UnitCellCoord(1, 0, 3, 1)},
      {UnitCellCoord(1, 1, -3, 0)},
      {UnitCellCoord(1, 1, 0, -3)},
      {UnitCellCoord(1, 1, 0, 3)},
      {UnitCellCoord(1, 1, 3, 0)},
      {UnitCellCoord(1, 3, -1, 0)},
      {UnitCellCoord(1, 3, 0, -1)},
      {UnitCellCoord(1, 3, 0, 1)},
      {UnitCellCoord(1, 3, 1, 0)}
    };

    m_orbit_neighborhood[12] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -3, -1, -1)},
      {UnitCellCoord(1, -3, -1, 1)},
      {UnitCellCoord(1, -3, 1, -1)},
      {UnitCellCoord(1, -3, 1, 1)},
      {UnitCellCoord(1, -1, -3, -1)},
      {UnitCellCoord(1, -1, -3, 1)},
      {UnitCellCoord(1, -1, -1, -3)},
      {UnitCellCoord(1, -1, -1, 3)},
      {UnitCellCoord(1, -1, 1, -3)},
      {UnitCellCoord(1, -1, 1, 3)},
      {UnitCellCoord(1, -1, 3, -1)},
      {UnitCellCoord(1, -1, 3, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 1, -3, -1)},
      {UnitCellCoord(1, 1, -3, 1)},
      {UnitCellCoord(1, 1, -1, -3)},
      {UnitCellCoord(1, 1, -1, 3)},
      {UnitCellCoord(1, 1, 1, -3)},
      {UnitCellCoord(1, 1, 1, 3)},
      {UnitCellCoord(1, 1, 3, -1)},
      {UnitCellCoord(1, 1, 3, 1)},
      {UnitCellCoord(1, 3, -1, -1)},
      {UnitCellCoord(1, 3, -1, 1)},
      {UnitCellCoord(1, 3, 1, -1)},
      {UnitCellCoord(1, 3, 1, 1)}
    };

    m_orbit_neighborhood[13] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -2, -2)},
      {UnitCellCoord(1, -2, -2, 2)},
      {UnitCellCoord(1, -2, 2, -2)},
      {UnitCellCoord(1, -2, 2, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 2, -2, -2)},
      {UnitCellCoord(1, 2, -2, 2)},
      {UnitCellCoord(1, 2, 2, -2)},
      {UnitCellCoord(1, 2, 2, 2)}
    };

    m_orbit_neighborhood[14] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, 0)}
    };

    m_orbit_neighborhood[15] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, 0)}
    };

    m_orbit_neighborhood[16] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 1, 1)}
    };

    m_orbit_neighborhood[17] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 2, 0, 0)}
    };

    m_orbit_neighborhood[18] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 2, 0, 0)}
    };

    m_orbit_neighborhood[19] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 1)},
      {UnitCellCoord(1, 2, 0, 0)}
    };

    m_orbit_neighborhood[20] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 0)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, 0)}
    };

    m_orbit_neighborhood[21] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, 0)}
    };

    m_orbit_neighborhood[22] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, 0)}
    };

    m_orbit_neighborhood[23] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 1, 1)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, 0)}
    };

    m_orbit_neighborhood[24] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, 0)}
    };

    m_orbit_neighborhood[25] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 0)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, 0)}
    };

    m_orbit_neighborhood[26] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 0)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[27] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 1)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[28] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[29] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 0, 0)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[30] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[31] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 1)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 0)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[32] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 0, 0)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[33] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 0)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

    m_orbit_neighborhood[34] = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 1)}
    };

  }

  La_MnNi_O3_cubic_Clexulator::~La_MnNi_O3_cubic_Clexulator(){
    //nothing here for now
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  void La_MnNi_O3_cubic_Clexulator::calc_global_corr_contribution(double *corr_begin) const {
    for(size_type i=0; i<corr_size(); i++){
      *(corr_begin+i) = (this->*m_orbit_func_list[i])();
    }
  }

  /// \brief Calculate contribution to select global correlations from one unit cell
  void La_MnNi_O3_cubic_Clexulator::calc_restricted_global_corr_contribution(double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {
    for(; ind_list_begin<ind_list_end; ind_list_begin++){
      *(corr_begin+*ind_list_begin) = (this->*m_orbit_func_list[*ind_list_begin])();
    }
  }

  /// \brief Calculate point correlations about basis site 'b_index'
  void La_MnNi_O3_cubic_Clexulator::calc_point_corr(int b_index, double *corr_begin) const {
    for(size_type i=0; i<corr_size(); i++){
      *(corr_begin+i) = (this->*m_flower_func_lists[b_index][i])();
    }
  }

  /// \brief Calculate select point correlations about basis site 'b_index'
  void La_MnNi_O3_cubic_Clexulator::calc_restricted_point_corr(int b_index, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {
    for(; ind_list_begin<ind_list_end; ind_list_begin++){
      *(corr_begin+*ind_list_begin) = (this->*m_flower_func_lists[b_index][*ind_list_begin])();
    }
  }

  /// \brief Calculate the change in point correlations due to changing an occupant
  void La_MnNi_O3_cubic_Clexulator::calc_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin) const {
    for(size_type i=0; i<corr_size(); i++){
      *(corr_begin+i) = (this->*m_delta_func_lists[b_index][i])(occ_i, occ_f);
    }
  }

  /// \brief Calculate the change in select point correlations due to changing an occupant
  void La_MnNi_O3_cubic_Clexulator::calc_restricted_delta_point_corr(int b_index, int occ_i, int occ_f, double *corr_begin, size_type const* ind_list_begin, size_type const* ind_list_end) const {
    for(; ind_list_begin<ind_list_end; ind_list_begin++){
      *(corr_begin+*ind_list_begin) = (this->*m_delta_func_lists[b_index][*ind_list_begin])(occ_i, occ_f);
    }
  }

  // Basis functions for empty cluster:
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_0_0_0() const{
    return (1);
  }

  /**** Basis functions for orbit 1, 0****
#Points: 1
MaxLength: 0  MinLength: 0
               0.5000000    0.5000000    0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_1_0_0() const{
    return (occ_func_1_0(0));
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_1_0_0() const{
    return (occ_func_1_0(0));
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_1_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i]);
  }

  /**** Basis functions for orbit 2, 0****
#Points: 2
MaxLength: 3.9401447  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000    0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_0_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(0)*occ_func_1_0(5)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_0_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(2)*occ_func_1_0(0)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(3)) + (occ_func_1_0(4)) + (occ_func_1_0(1)) + (occ_func_1_0(6)) + (occ_func_1_0(5)) + (occ_func_1_0(2)))/3.0;
  }

  /**** Basis functions for orbit 2, 1****
#Points: 2
MaxLength: 5.5722061  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_1_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(15)))/6.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_1_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(10)*occ_func_1_0(0)))/6.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(11)) + (occ_func_1_0(14)) + (occ_func_1_0(16)) + (occ_func_1_0(9)) + (occ_func_1_0(7)) + (occ_func_1_0(18)) + (occ_func_1_0(13)) + (occ_func_1_0(12)) + (occ_func_1_0(8)) + (occ_func_1_0(17)) + (occ_func_1_0(15)) + (occ_func_1_0(10)))/6.0;
  }

  /**** Basis functions for orbit 2, 2****
#Points: 2
MaxLength: 6.8245308  MinLength: 6.8245308
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_2_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(21)))/4.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_2_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(24)*occ_func_1_0(0)))/4.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(19)) + (occ_func_1_0(26)) + (occ_func_1_0(23)) + (occ_func_1_0(22)) + (occ_func_1_0(20)) + (occ_func_1_0(25)) + (occ_func_1_0(21)) + (occ_func_1_0(24)))/4.0;
  }

  /**** Basis functions for orbit 2, 3****
#Points: 2
MaxLength: 7.8802894  MinLength: 7.8802894
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000    0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_3_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(31)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_3_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)) + (occ_func_1_0(30)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(28)*occ_func_1_0(0)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(29)) + (occ_func_1_0(30)) + (occ_func_1_0(27)) + (occ_func_1_0(32)) + (occ_func_1_0(31)) + (occ_func_1_0(28)))/3.0;
  }

  /**** Basis functions for orbit 2, 4****
#Points: 2
MaxLength: 8.8104314  MinLength: 8.8104314
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_4_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(48)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_4_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)) + (occ_func_1_0(46)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)) + (occ_func_1_0(39)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)) + (occ_func_1_0(56)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)) + (occ_func_1_0(42)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(38)) + (occ_func_1_0(51)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)) + (occ_func_1_0(36)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(54)) + (occ_func_1_0(35)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)) + (occ_func_1_0(37)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(40)) + (occ_func_1_0(49)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(34)) + (occ_func_1_0(55)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(45)) + (occ_func_1_0(44)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(48)) + (occ_func_1_0(41)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_4_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(43)) + (occ_func_1_0(46)) + (occ_func_1_0(50)) + (occ_func_1_0(39)) + (occ_func_1_0(33)) + (occ_func_1_0(56)) + (occ_func_1_0(47)) + (occ_func_1_0(42)) + (occ_func_1_0(38)) + (occ_func_1_0(51)) + (occ_func_1_0(53)) + (occ_func_1_0(36)) + (occ_func_1_0(54)) + (occ_func_1_0(35)) + (occ_func_1_0(52)) + (occ_func_1_0(37)) + (occ_func_1_0(40)) + (occ_func_1_0(49)) + (occ_func_1_0(34)) + (occ_func_1_0(55)) + (occ_func_1_0(45)) + (occ_func_1_0(44)) + (occ_func_1_0(48)) + (occ_func_1_0(41)))/12.0;
  }

  /**** Basis functions for orbit 2, 5****
#Points: 2
MaxLength: 9.6513440  MinLength: 9.6513440
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_5_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)) + (occ_func_1_0(0)*occ_func_1_0(71)) + (occ_func_1_0(0)*occ_func_1_0(58)) + (occ_func_1_0(0)*occ_func_1_0(67)) + (occ_func_1_0(0)*occ_func_1_0(65)) + (occ_func_1_0(0)*occ_func_1_0(77)) + (occ_func_1_0(0)*occ_func_1_0(62)) + (occ_func_1_0(0)*occ_func_1_0(68)) + (occ_func_1_0(0)*occ_func_1_0(57)) + (occ_func_1_0(0)*occ_func_1_0(78)) + (occ_func_1_0(0)*occ_func_1_0(61)) + (occ_func_1_0(0)*occ_func_1_0(73)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_5_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)) + (occ_func_1_0(74)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)) + (occ_func_1_0(66)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)) + (occ_func_1_0(79)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)) + (occ_func_1_0(70)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)) + (occ_func_1_0(72)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)) + (occ_func_1_0(60)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)) + (occ_func_1_0(75)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)) + (occ_func_1_0(69)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)) + (occ_func_1_0(80)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)) + (occ_func_1_0(59)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)) + (occ_func_1_0(76)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)) + (occ_func_1_0(64)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_5_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)) + (occ_func_1_0(74)) + (occ_func_1_0(71)) + (occ_func_1_0(66)) + (occ_func_1_0(58)) + (occ_func_1_0(79)) + (occ_func_1_0(67)) + (occ_func_1_0(70)) + (occ_func_1_0(65)) + (occ_func_1_0(72)) + (occ_func_1_0(77)) + (occ_func_1_0(60)) + (occ_func_1_0(62)) + (occ_func_1_0(75)) + (occ_func_1_0(68)) + (occ_func_1_0(69)) + (occ_func_1_0(57)) + (occ_func_1_0(80)) + (occ_func_1_0(78)) + (occ_func_1_0(59)) + (occ_func_1_0(61)) + (occ_func_1_0(76)) + (occ_func_1_0(73)) + (occ_func_1_0(64)))/12.0;
  }

  /**** Basis functions for orbit 2, 6****
#Points: 2
MaxLength: 11.1444122  MinLength: 11.1444122
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -1.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_6_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(85)) + (occ_func_1_0(0)*occ_func_1_0(90)) + (occ_func_1_0(0)*occ_func_1_0(81)) + (occ_func_1_0(0)*occ_func_1_0(87)) + (occ_func_1_0(0)*occ_func_1_0(82)) + (occ_func_1_0(0)*occ_func_1_0(89)))/6.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_6_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(85)) + (occ_func_1_0(88)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(90)) + (occ_func_1_0(83)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(81)) + (occ_func_1_0(92)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(87)) + (occ_func_1_0(86)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(82)) + (occ_func_1_0(91)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(89)) + (occ_func_1_0(84)*occ_func_1_0(0)))/6.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_6_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(85)) + (occ_func_1_0(88)) + (occ_func_1_0(90)) + (occ_func_1_0(83)) + (occ_func_1_0(81)) + (occ_func_1_0(92)) + (occ_func_1_0(87)) + (occ_func_1_0(86)) + (occ_func_1_0(82)) + (occ_func_1_0(91)) + (occ_func_1_0(89)) + (occ_func_1_0(84)))/6.0;
  }

  /**** Basis functions for orbit 2, 7****
#Points: 2
MaxLength: 11.8204341  MinLength: 11.8204341
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -1.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_7_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(102)) + (occ_func_1_0(0)*occ_func_1_0(116)) + (occ_func_1_0(0)*occ_func_1_0(95)) + (occ_func_1_0(0)*occ_func_1_0(104)) + (occ_func_1_0(0)*occ_func_1_0(98)) + (occ_func_1_0(0)*occ_func_1_0(114)) + (occ_func_1_0(0)*occ_func_1_0(103)) + (occ_func_1_0(0)*occ_func_1_0(118)) + (occ_func_1_0(0)*occ_func_1_0(96)) + (occ_func_1_0(0)*occ_func_1_0(115)) + (occ_func_1_0(0)*occ_func_1_0(94)) + (occ_func_1_0(0)*occ_func_1_0(105)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_7_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(102)) + (occ_func_1_0(113)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(116)) + (occ_func_1_0(99)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(95)) + (occ_func_1_0(120)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(104)) + (occ_func_1_0(111)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(98)) + (occ_func_1_0(117)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(114)) + (occ_func_1_0(101)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(103)) + (occ_func_1_0(112)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(118)) + (occ_func_1_0(97)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(96)) + (occ_func_1_0(119)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(115)) + (occ_func_1_0(100)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(94)) + (occ_func_1_0(121)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(105)) + (occ_func_1_0(110)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_7_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(102)) + (occ_func_1_0(113)) + (occ_func_1_0(116)) + (occ_func_1_0(99)) + (occ_func_1_0(95)) + (occ_func_1_0(120)) + (occ_func_1_0(104)) + (occ_func_1_0(111)) + (occ_func_1_0(98)) + (occ_func_1_0(117)) + (occ_func_1_0(114)) + (occ_func_1_0(101)) + (occ_func_1_0(103)) + (occ_func_1_0(112)) + (occ_func_1_0(118)) + (occ_func_1_0(97)) + (occ_func_1_0(96)) + (occ_func_1_0(119)) + (occ_func_1_0(115)) + (occ_func_1_0(100)) + (occ_func_1_0(94)) + (occ_func_1_0(121)) + (occ_func_1_0(105)) + (occ_func_1_0(110)))/12.0;
  }

  /**** Basis functions for orbit 2, 8****
#Points: 2
MaxLength: 11.8204341  MinLength: 11.8204341
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000    0.5000000   -2.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_8_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(107)) + (occ_func_1_0(0)*occ_func_1_0(93)) + (occ_func_1_0(0)*occ_func_1_0(109)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_8_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(107)) + (occ_func_1_0(108)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(93)) + (occ_func_1_0(122)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(109)) + (occ_func_1_0(106)*occ_func_1_0(0)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_8_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(107)) + (occ_func_1_0(108)) + (occ_func_1_0(93)) + (occ_func_1_0(122)) + (occ_func_1_0(109)) + (occ_func_1_0(106)))/3.0;
  }

  /**** Basis functions for orbit 2, 9****
#Points: 2
MaxLength: 12.4598316  MinLength: 12.4598316
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -2.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_9_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(133)) + (occ_func_1_0(0)*occ_func_1_0(140)) + (occ_func_1_0(0)*occ_func_1_0(123)) + (occ_func_1_0(0)*occ_func_1_0(137)) + (occ_func_1_0(0)*occ_func_1_0(128)) + (occ_func_1_0(0)*occ_func_1_0(143)) + (occ_func_1_0(0)*occ_func_1_0(144)) + (occ_func_1_0(0)*occ_func_1_0(142)) + (occ_func_1_0(0)*occ_func_1_0(130)) + (occ_func_1_0(0)*occ_func_1_0(124)) + (occ_func_1_0(0)*occ_func_1_0(135)) + (occ_func_1_0(0)*occ_func_1_0(138)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_9_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(133)) + (occ_func_1_0(136)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(140)) + (occ_func_1_0(129)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(123)) + (occ_func_1_0(146)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(137)) + (occ_func_1_0(132)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(128)) + (occ_func_1_0(141)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(143)) + (occ_func_1_0(126)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(144)) + (occ_func_1_0(125)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(142)) + (occ_func_1_0(127)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(130)) + (occ_func_1_0(139)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(124)) + (occ_func_1_0(145)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(135)) + (occ_func_1_0(134)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(138)) + (occ_func_1_0(131)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_9_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(133)) + (occ_func_1_0(136)) + (occ_func_1_0(140)) + (occ_func_1_0(129)) + (occ_func_1_0(123)) + (occ_func_1_0(146)) + (occ_func_1_0(137)) + (occ_func_1_0(132)) + (occ_func_1_0(128)) + (occ_func_1_0(141)) + (occ_func_1_0(143)) + (occ_func_1_0(126)) + (occ_func_1_0(144)) + (occ_func_1_0(125)) + (occ_func_1_0(142)) + (occ_func_1_0(127)) + (occ_func_1_0(130)) + (occ_func_1_0(139)) + (occ_func_1_0(124)) + (occ_func_1_0(145)) + (occ_func_1_0(135)) + (occ_func_1_0(134)) + (occ_func_1_0(138)) + (occ_func_1_0(131)))/12.0;
  }

  /**** Basis functions for orbit 2, 10****
#Points: 2
MaxLength: 13.0679816  MinLength: 13.0679816
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -2.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_10_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(153)) + (occ_func_1_0(0)*occ_func_1_0(161)) + (occ_func_1_0(0)*occ_func_1_0(148)) + (occ_func_1_0(0)*occ_func_1_0(157)) + (occ_func_1_0(0)*occ_func_1_0(155)) + (occ_func_1_0(0)*occ_func_1_0(167)) + (occ_func_1_0(0)*occ_func_1_0(152)) + (occ_func_1_0(0)*occ_func_1_0(158)) + (occ_func_1_0(0)*occ_func_1_0(147)) + (occ_func_1_0(0)*occ_func_1_0(168)) + (occ_func_1_0(0)*occ_func_1_0(151)) + (occ_func_1_0(0)*occ_func_1_0(163)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_10_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(153)) + (occ_func_1_0(164)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(161)) + (occ_func_1_0(156)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(148)) + (occ_func_1_0(169)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(157)) + (occ_func_1_0(160)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(155)) + (occ_func_1_0(162)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(167)) + (occ_func_1_0(150)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(152)) + (occ_func_1_0(165)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(158)) + (occ_func_1_0(159)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(147)) + (occ_func_1_0(170)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(168)) + (occ_func_1_0(149)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(151)) + (occ_func_1_0(166)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(163)) + (occ_func_1_0(154)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_10_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(153)) + (occ_func_1_0(164)) + (occ_func_1_0(161)) + (occ_func_1_0(156)) + (occ_func_1_0(148)) + (occ_func_1_0(169)) + (occ_func_1_0(157)) + (occ_func_1_0(160)) + (occ_func_1_0(155)) + (occ_func_1_0(162)) + (occ_func_1_0(167)) + (occ_func_1_0(150)) + (occ_func_1_0(152)) + (occ_func_1_0(165)) + (occ_func_1_0(158)) + (occ_func_1_0(159)) + (occ_func_1_0(147)) + (occ_func_1_0(170)) + (occ_func_1_0(168)) + (occ_func_1_0(149)) + (occ_func_1_0(151)) + (occ_func_1_0(166)) + (occ_func_1_0(163)) + (occ_func_1_0(154)))/12.0;
  }

  /**** Basis functions for orbit 2, 11****
#Points: 2
MaxLength: 13.6490616  MinLength: 13.6490616
               0.5000000    0.5000000    0.5000000 Mn Ni
              -1.5000000   -1.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_11_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(171)) + (occ_func_1_0(0)*occ_func_1_0(175)) + (occ_func_1_0(0)*occ_func_1_0(172)) + (occ_func_1_0(0)*occ_func_1_0(173)))/4.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_11_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(171)) + (occ_func_1_0(178)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(175)) + (occ_func_1_0(174)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(172)) + (occ_func_1_0(177)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(173)) + (occ_func_1_0(176)*occ_func_1_0(0)))/4.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_11_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(171)) + (occ_func_1_0(178)) + (occ_func_1_0(175)) + (occ_func_1_0(174)) + (occ_func_1_0(172)) + (occ_func_1_0(177)) + (occ_func_1_0(173)) + (occ_func_1_0(176)))/4.0;
  }

  /**** Basis functions for orbit 3, 0****
#Points: 3
MaxLength: 5.5722061  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -0.5000000 Mn Ni
               0.5000000    0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_0_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(11)*occ_func_1_0(3)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(3)) + (occ_func_1_0(0)*occ_func_1_0(7)*occ_func_1_0(1)) + (occ_func_1_0(0)*occ_func_1_0(13)*occ_func_1_0(5)) + (occ_func_1_0(0)*occ_func_1_0(8)*occ_func_1_0(3)) + (occ_func_1_0(0)*occ_func_1_0(15)*occ_func_1_0(6)) + (occ_func_1_0(0)*occ_func_1_0(12)*occ_func_1_0(2)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(6)) + (occ_func_1_0(0)*occ_func_1_0(8)*occ_func_1_0(1)) + (occ_func_1_0(0)*occ_func_1_0(15)*occ_func_1_0(2)) + (occ_func_1_0(0)*occ_func_1_0(7)*occ_func_1_0(2)) + (occ_func_1_0(0)*occ_func_1_0(14)*occ_func_1_0(4)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_0_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(11)*occ_func_1_0(3)) + (occ_func_1_0(14)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(4)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(3)) + (occ_func_1_0(9)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(4)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(7)*occ_func_1_0(1)) + (occ_func_1_0(18)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(6)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(13)*occ_func_1_0(5)) + (occ_func_1_0(12)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(2)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(8)*occ_func_1_0(3)) + (occ_func_1_0(17)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(4)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(15)*occ_func_1_0(6)) + (occ_func_1_0(10)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(1)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(12)*occ_func_1_0(2)) + (occ_func_1_0(13)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(5)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(6)) + (occ_func_1_0(9)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(1)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(8)*occ_func_1_0(1)) + (occ_func_1_0(17)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(6)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(15)*occ_func_1_0(2)) + (occ_func_1_0(10)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(5)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(7)*occ_func_1_0(2)) + (occ_func_1_0(18)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(5)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(14)*occ_func_1_0(4)) + (occ_func_1_0(11)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(3)*occ_func_1_0(5)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_0_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(11)*occ_func_1_0(3)) + (occ_func_1_0(14)*occ_func_1_0(5)) + (occ_func_1_0(4)*occ_func_1_0(2)) + (occ_func_1_0(16)*occ_func_1_0(3)) + (occ_func_1_0(9)*occ_func_1_0(1)) + (occ_func_1_0(4)*occ_func_1_0(6)) + (occ_func_1_0(7)*occ_func_1_0(1)) + (occ_func_1_0(18)*occ_func_1_0(5)) + (occ_func_1_0(6)*occ_func_1_0(2)) + (occ_func_1_0(13)*occ_func_1_0(5)) + (occ_func_1_0(12)*occ_func_1_0(4)) + (occ_func_1_0(2)*occ_func_1_0(3)) + (occ_func_1_0(8)*occ_func_1_0(3)) + (occ_func_1_0(17)*occ_func_1_0(6)) + (occ_func_1_0(4)*occ_func_1_0(1)) + (occ_func_1_0(15)*occ_func_1_0(6)) + (occ_func_1_0(10)*occ_func_1_0(5)) + (occ_func_1_0(1)*occ_func_1_0(2)) + (occ_func_1_0(12)*occ_func_1_0(2)) + (occ_func_1_0(13)*occ_func_1_0(3)) + (occ_func_1_0(5)*occ_func_1_0(4)) + (occ_func_1_0(16)*occ_func_1_0(6)) + (occ_func_1_0(9)*occ_func_1_0(4)) + (occ_func_1_0(1)*occ_func_1_0(3)) + (occ_func_1_0(8)*occ_func_1_0(1)) + (occ_func_1_0(17)*occ_func_1_0(4)) + (occ_func_1_0(6)*occ_func_1_0(3)) + (occ_func_1_0(15)*occ_func_1_0(2)) + (occ_func_1_0(10)*occ_func_1_0(1)) + (occ_func_1_0(5)*occ_func_1_0(6)) + (occ_func_1_0(7)*occ_func_1_0(2)) + (occ_func_1_0(18)*occ_func_1_0(6)) + (occ_func_1_0(5)*occ_func_1_0(1)) + (occ_func_1_0(14)*occ_func_1_0(4)) + (occ_func_1_0(11)*occ_func_1_0(2)) + (occ_func_1_0(3)*occ_func_1_0(5)))/12.0;
  }

  /**** Basis functions for orbit 3, 1****
#Points: 3
MaxLength: 5.5722061  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -0.5000000 Mn Ni
              -0.5000000    0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_1_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(11)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(13)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(8)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(18)) + (occ_func_1_0(0)*occ_func_1_0(18)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(15)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(11)*occ_func_1_0(15)))/8.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_1_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(11)*occ_func_1_0(8)) + (occ_func_1_0(14)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(17)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(11)) + (occ_func_1_0(9)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(14)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(13)*occ_func_1_0(10)) + (occ_func_1_0(12)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(15)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(8)*occ_func_1_0(13)) + (occ_func_1_0(17)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(12)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(16)*occ_func_1_0(18)) + (occ_func_1_0(9)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(7)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(18)*occ_func_1_0(13)) + (occ_func_1_0(7)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(12)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(15)*occ_func_1_0(12)) + (occ_func_1_0(10)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(13)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(11)*occ_func_1_0(15)) + (occ_func_1_0(14)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(10)*occ_func_1_0(8)*occ_func_1_0(0)))/8.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_1_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(11)*occ_func_1_0(8)) + (occ_func_1_0(14)*occ_func_1_0(10)) + (occ_func_1_0(17)*occ_func_1_0(15)) + (occ_func_1_0(16)*occ_func_1_0(11)) + (occ_func_1_0(9)*occ_func_1_0(7)) + (occ_func_1_0(14)*occ_func_1_0(18)) + (occ_func_1_0(13)*occ_func_1_0(10)) + (occ_func_1_0(12)*occ_func_1_0(9)) + (occ_func_1_0(15)*occ_func_1_0(16)) + (occ_func_1_0(8)*occ_func_1_0(13)) + (occ_func_1_0(17)*occ_func_1_0(18)) + (occ_func_1_0(12)*occ_func_1_0(7)) + (occ_func_1_0(16)*occ_func_1_0(18)) + (occ_func_1_0(9)*occ_func_1_0(14)) + (occ_func_1_0(7)*occ_func_1_0(11)) + (occ_func_1_0(18)*occ_func_1_0(13)) + (occ_func_1_0(7)*occ_func_1_0(8)) + (occ_func_1_0(12)*occ_func_1_0(17)) + (occ_func_1_0(15)*occ_func_1_0(12)) + (occ_func_1_0(10)*occ_func_1_0(9)) + (occ_func_1_0(13)*occ_func_1_0(16)) + (occ_func_1_0(11)*occ_func_1_0(15)) + (occ_func_1_0(14)*occ_func_1_0(17)) + (occ_func_1_0(10)*occ_func_1_0(8)))/8.0;
  }

  /**** Basis functions for orbit 3, 2****
#Points: 3
MaxLength: 6.8245308  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -0.5000000 Mn Ni
               0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_2_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(19)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(23)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(20)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(21)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(21)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(23)*occ_func_1_0(15)) + (occ_func_1_0(0)*occ_func_1_0(20)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(25)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(25)*occ_func_1_0(18)) + (occ_func_1_0(0)*occ_func_1_0(22)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(19)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(24)*occ_func_1_0(15)) + (occ_func_1_0(0)*occ_func_1_0(22)*occ_func_1_0(9)) + (occ_func_1_0(0)*occ_func_1_0(24)*occ_func_1_0(17)) + (occ_func_1_0(0)*occ_func_1_0(19)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(25)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(24)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(22)*occ_func_1_0(14)) + (occ_func_1_0(0)*occ_func_1_0(26)*occ_func_1_0(14)) + (occ_func_1_0(0)*occ_func_1_0(23)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(26)*occ_func_1_0(17)) + (occ_func_1_0(0)*occ_func_1_0(26)*occ_func_1_0(18)) + (occ_func_1_0(0)*occ_func_1_0(21)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(20)*occ_func_1_0(9)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_2_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(19)*occ_func_1_0(11)) + (occ_func_1_0(26)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(14)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(23)*occ_func_1_0(16)) + (occ_func_1_0(22)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(9)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(20)*occ_func_1_0(7)) + (occ_func_1_0(25)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(18)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(21)*occ_func_1_0(13)) + (occ_func_1_0(24)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(12)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(21)*occ_func_1_0(8)) + (occ_func_1_0(24)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(17)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(23)*occ_func_1_0(15)) + (occ_func_1_0(22)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(10)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(20)*occ_func_1_0(12)) + (occ_func_1_0(25)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(13)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(25)*occ_func_1_0(16)) + (occ_func_1_0(20)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(9)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(25)*occ_func_1_0(18)) + (occ_func_1_0(20)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(7)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(22)*occ_func_1_0(10)) + (occ_func_1_0(23)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(15)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(19)*occ_func_1_0(8)) + (occ_func_1_0(26)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(17)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(24)*occ_func_1_0(15)) + (occ_func_1_0(21)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(10)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(22)*occ_func_1_0(9)) + (occ_func_1_0(23)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(16)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(24)*occ_func_1_0(17)) + (occ_func_1_0(21)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(8)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(19)*occ_func_1_0(7)) + (occ_func_1_0(26)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(18)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(25)*occ_func_1_0(13)) + (occ_func_1_0(20)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(12)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(24)*occ_func_1_0(12)) + (occ_func_1_0(21)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(13)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(22)*occ_func_1_0(14)) + (occ_func_1_0(23)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(11)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(26)*occ_func_1_0(14)) + (occ_func_1_0(19)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(11)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(23)*occ_func_1_0(11)) + (occ_func_1_0(22)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(14)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(26)*occ_func_1_0(17)) + (occ_func_1_0(19)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(8)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(26)*occ_func_1_0(18)) + (occ_func_1_0(19)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(7)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(21)*occ_func_1_0(10)) + (occ_func_1_0(24)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(15)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(20)*occ_func_1_0(9)) + (occ_func_1_0(25)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(16)*occ_func_1_0(2)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_2_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(19)*occ_func_1_0(11)) + (occ_func_1_0(26)*occ_func_1_0(6)) + (occ_func_1_0(14)*occ_func_1_0(1)) + (occ_func_1_0(23)*occ_func_1_0(16)) + (occ_func_1_0(22)*occ_func_1_0(5)) + (occ_func_1_0(9)*occ_func_1_0(2)) + (occ_func_1_0(20)*occ_func_1_0(7)) + (occ_func_1_0(25)*occ_func_1_0(3)) + (occ_func_1_0(18)*occ_func_1_0(4)) + (occ_func_1_0(21)*occ_func_1_0(13)) + (occ_func_1_0(24)*occ_func_1_0(6)) + (occ_func_1_0(12)*occ_func_1_0(1)) + (occ_func_1_0(21)*occ_func_1_0(8)) + (occ_func_1_0(24)*occ_func_1_0(2)) + (occ_func_1_0(17)*occ_func_1_0(5)) + (occ_func_1_0(23)*occ_func_1_0(15)) + (occ_func_1_0(22)*occ_func_1_0(4)) + (occ_func_1_0(10)*occ_func_1_0(3)) + (occ_func_1_0(20)*occ_func_1_0(12)) + (occ_func_1_0(25)*occ_func_1_0(6)) + (occ_func_1_0(13)*occ_func_1_0(1)) + (occ_func_1_0(25)*occ_func_1_0(16)) + (occ_func_1_0(20)*occ_func_1_0(2)) + (occ_func_1_0(9)*occ_func_1_0(5)) + (occ_func_1_0(25)*occ_func_1_0(18)) + (occ_func_1_0(20)*occ_func_1_0(4)) + (occ_func_1_0(7)*occ_func_1_0(3)) + (occ_func_1_0(22)*occ_func_1_0(10)) + (occ_func_1_0(23)*occ_func_1_0(3)) + (occ_func_1_0(15)*occ_func_1_0(4)) + (occ_func_1_0(19)*occ_func_1_0(8)) + (occ_func_1_0(26)*occ_func_1_0(5)) + (occ_func_1_0(17)*occ_func_1_0(2)) + (occ_func_1_0(24)*occ_func_1_0(15)) + (occ_func_1_0(21)*occ_func_1_0(3)) + (occ_func_1_0(10)*occ_func_1_0(4)) + (occ_func_1_0(22)*occ_func_1_0(9)) + (occ_func_1_0(23)*occ_func_1_0(2)) + (occ_func_1_0(16)*occ_func_1_0(5)) + (occ_func_1_0(24)*occ_func_1_0(17)) + (occ_func_1_0(21)*occ_func_1_0(5)) + (occ_func_1_0(8)*occ_func_1_0(2)) + (occ_func_1_0(19)*occ_func_1_0(7)) + (occ_func_1_0(26)*occ_func_1_0(4)) + (occ_func_1_0(18)*occ_func_1_0(3)) + (occ_func_1_0(25)*occ_func_1_0(13)) + (occ_func_1_0(20)*occ_func_1_0(1)) + (occ_func_1_0(12)*occ_func_1_0(6)) + (occ_func_1_0(24)*occ_func_1_0(12)) + (occ_func_1_0(21)*occ_func_1_0(1)) + (occ_func_1_0(13)*occ_func_1_0(6)) + (occ_func_1_0(22)*occ_func_1_0(14)) + (occ_func_1_0(23)*occ_func_1_0(6)) + (occ_func_1_0(11)*occ_func_1_0(1)) + (occ_func_1_0(26)*occ_func_1_0(14)) + (occ_func_1_0(19)*occ_func_1_0(1)) + (occ_func_1_0(11)*occ_func_1_0(6)) + (occ_func_1_0(23)*occ_func_1_0(11)) + (occ_func_1_0(22)*occ_func_1_0(1)) + (occ_func_1_0(14)*occ_func_1_0(6)) + (occ_func_1_0(26)*occ_func_1_0(17)) + (occ_func_1_0(19)*occ_func_1_0(2)) + (occ_func_1_0(8)*occ_func_1_0(5)) + (occ_func_1_0(26)*occ_func_1_0(18)) + (occ_func_1_0(19)*occ_func_1_0(3)) + (occ_func_1_0(7)*occ_func_1_0(4)) + (occ_func_1_0(21)*occ_func_1_0(10)) + (occ_func_1_0(24)*occ_func_1_0(4)) + (occ_func_1_0(15)*occ_func_1_0(3)) + (occ_func_1_0(20)*occ_func_1_0(9)) + (occ_func_1_0(25)*occ_func_1_0(5)) + (occ_func_1_0(16)*occ_func_1_0(2)))/24.0;
  }

  /**** Basis functions for orbit 3, 3****
#Points: 3
MaxLength: 7.8802894  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000    0.5000000   -1.5000000 Mn Ni
               0.5000000    0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_3_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(3)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(1)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(5)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_3_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(3)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(4)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(1)) + (occ_func_1_0(32)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(6)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(5)) + (occ_func_1_0(28)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(2)*occ_func_1_0(5)*occ_func_1_0(0)))/3.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_3_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(29)*occ_func_1_0(3)) + (occ_func_1_0(30)*occ_func_1_0(4)) + (occ_func_1_0(4)*occ_func_1_0(3)) + (occ_func_1_0(27)*occ_func_1_0(1)) + (occ_func_1_0(32)*occ_func_1_0(6)) + (occ_func_1_0(6)*occ_func_1_0(1)) + (occ_func_1_0(31)*occ_func_1_0(5)) + (occ_func_1_0(28)*occ_func_1_0(2)) + (occ_func_1_0(2)*occ_func_1_0(5)))/3.0;
  }

  /**** Basis functions for orbit 3, 4****
#Points: 3
MaxLength: 7.8802894  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000    0.5000000   -1.5000000 Mn Ni
               0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_4_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(28)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(18)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(9)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(18)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_4_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(11)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(14)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(16)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(9)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(7)) + (occ_func_1_0(32)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(18)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(13)) + (occ_func_1_0(28)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(12)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(8)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(17)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(28)*occ_func_1_0(12)) + (occ_func_1_0(31)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(13)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(16)) + (occ_func_1_0(27)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(9)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(18)) + (occ_func_1_0(28)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(7)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(10)) + (occ_func_1_0(28)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(15)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(9)) + (occ_func_1_0(32)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(16)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(13)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(12)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(18)) + (occ_func_1_0(27)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(7)*occ_func_1_0(15)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_4_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(29)*occ_func_1_0(11)) + (occ_func_1_0(30)*occ_func_1_0(12)) + (occ_func_1_0(14)*occ_func_1_0(13)) + (occ_func_1_0(29)*occ_func_1_0(16)) + (occ_func_1_0(30)*occ_func_1_0(17)) + (occ_func_1_0(9)*occ_func_1_0(8)) + (occ_func_1_0(27)*occ_func_1_0(7)) + (occ_func_1_0(32)*occ_func_1_0(15)) + (occ_func_1_0(18)*occ_func_1_0(10)) + (occ_func_1_0(31)*occ_func_1_0(13)) + (occ_func_1_0(28)*occ_func_1_0(11)) + (occ_func_1_0(12)*occ_func_1_0(14)) + (occ_func_1_0(29)*occ_func_1_0(8)) + (occ_func_1_0(30)*occ_func_1_0(9)) + (occ_func_1_0(17)*occ_func_1_0(16)) + (occ_func_1_0(28)*occ_func_1_0(12)) + (occ_func_1_0(31)*occ_func_1_0(14)) + (occ_func_1_0(13)*occ_func_1_0(11)) + (occ_func_1_0(32)*occ_func_1_0(16)) + (occ_func_1_0(27)*occ_func_1_0(8)) + (occ_func_1_0(9)*occ_func_1_0(17)) + (occ_func_1_0(31)*occ_func_1_0(18)) + (occ_func_1_0(28)*occ_func_1_0(15)) + (occ_func_1_0(7)*occ_func_1_0(10)) + (occ_func_1_0(31)*occ_func_1_0(10)) + (occ_func_1_0(28)*occ_func_1_0(7)) + (occ_func_1_0(15)*occ_func_1_0(18)) + (occ_func_1_0(27)*occ_func_1_0(9)) + (occ_func_1_0(32)*occ_func_1_0(17)) + (occ_func_1_0(16)*occ_func_1_0(8)) + (occ_func_1_0(29)*occ_func_1_0(13)) + (occ_func_1_0(30)*occ_func_1_0(14)) + (occ_func_1_0(12)*occ_func_1_0(11)) + (occ_func_1_0(32)*occ_func_1_0(18)) + (occ_func_1_0(27)*occ_func_1_0(10)) + (occ_func_1_0(7)*occ_func_1_0(15)))/12.0;
  }

  /**** Basis functions for orbit 3, 5****
#Points: 3
MaxLength: 7.8802894  MinLength: 6.8245308
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000    0.5000000   -1.5000000 Mn Ni
              -0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_5_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(28)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(28)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(25)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_5_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(19)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(26)*occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(23)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(22)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(20)) + (occ_func_1_0(32)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(25)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(21)) + (occ_func_1_0(28)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(24)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(21)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(24)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(23)) + (occ_func_1_0(27)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(22)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(28)*occ_func_1_0(20)) + (occ_func_1_0(31)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(25)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(32)*occ_func_1_0(25)) + (occ_func_1_0(27)*occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(20)*occ_func_1_0(24)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(31)*occ_func_1_0(25)) + (occ_func_1_0(28)*occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(20)*occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(28)*occ_func_1_0(24)) + (occ_func_1_0(31)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(21)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(27)*occ_func_1_0(22)) + (occ_func_1_0(32)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(23)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(29)*occ_func_1_0(25)) + (occ_func_1_0(30)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(20)*occ_func_1_0(19)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_5_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(29)*occ_func_1_0(19)) + (occ_func_1_0(30)*occ_func_1_0(20)) + (occ_func_1_0(26)*occ_func_1_0(25)) + (occ_func_1_0(29)*occ_func_1_0(23)) + (occ_func_1_0(30)*occ_func_1_0(24)) + (occ_func_1_0(22)*occ_func_1_0(21)) + (occ_func_1_0(27)*occ_func_1_0(20)) + (occ_func_1_0(32)*occ_func_1_0(24)) + (occ_func_1_0(25)*occ_func_1_0(21)) + (occ_func_1_0(31)*occ_func_1_0(21)) + (occ_func_1_0(28)*occ_func_1_0(19)) + (occ_func_1_0(24)*occ_func_1_0(26)) + (occ_func_1_0(29)*occ_func_1_0(21)) + (occ_func_1_0(30)*occ_func_1_0(22)) + (occ_func_1_0(24)*occ_func_1_0(23)) + (occ_func_1_0(32)*occ_func_1_0(23)) + (occ_func_1_0(27)*occ_func_1_0(19)) + (occ_func_1_0(22)*occ_func_1_0(26)) + (occ_func_1_0(28)*occ_func_1_0(20)) + (occ_func_1_0(31)*occ_func_1_0(22)) + (occ_func_1_0(25)*occ_func_1_0(23)) + (occ_func_1_0(32)*occ_func_1_0(25)) + (occ_func_1_0(27)*occ_func_1_0(21)) + (occ_func_1_0(20)*occ_func_1_0(24)) + (occ_func_1_0(31)*occ_func_1_0(25)) + (occ_func_1_0(28)*occ_func_1_0(23)) + (occ_func_1_0(20)*occ_func_1_0(22)) + (occ_func_1_0(28)*occ_func_1_0(24)) + (occ_func_1_0(31)*occ_func_1_0(26)) + (occ_func_1_0(21)*occ_func_1_0(19)) + (occ_func_1_0(27)*occ_func_1_0(22)) + (occ_func_1_0(32)*occ_func_1_0(26)) + (occ_func_1_0(23)*occ_func_1_0(19)) + (occ_func_1_0(29)*occ_func_1_0(25)) + (occ_func_1_0(30)*occ_func_1_0(26)) + (occ_func_1_0(20)*occ_func_1_0(19)))/12.0;
  }

  /**** Basis functions for orbit 3, 6****
#Points: 3
MaxLength: 8.8104314  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000    0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_6_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(30)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(30)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(30)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(30)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_6_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(29)) + (occ_func_1_0(46)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(30)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(29)) + (occ_func_1_0(39)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(30)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(27)) + (occ_func_1_0(56)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(32)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(31)) + (occ_func_1_0(42)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(28)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(29)) + (occ_func_1_0(51)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(30)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(32)) + (occ_func_1_0(36)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(27)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(28)) + (occ_func_1_0(47)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(31)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(32)) + (occ_func_1_0(35)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(27)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(31)) + (occ_func_1_0(37)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(28)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(31)) + (occ_func_1_0(49)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(28)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(27)) + (occ_func_1_0(55)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(32)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(28)) + (occ_func_1_0(40)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(31)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(27)) + (occ_func_1_0(54)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(32)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(32)) + (occ_func_1_0(34)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(27)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(28)) + (occ_func_1_0(52)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(31)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(29)) + (occ_func_1_0(44)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(30)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(30)) + (occ_func_1_0(45)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(29)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(30)) + (occ_func_1_0(43)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(29)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(31)) + (occ_func_1_0(41)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(28)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(28)) + (occ_func_1_0(48)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(31)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(30)) + (occ_func_1_0(38)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(29)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(32)) + (occ_func_1_0(33)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(27)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(27)) + (occ_func_1_0(53)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(32)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(30)) + (occ_func_1_0(50)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(29)*occ_func_1_0(1)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_6_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(43)*occ_func_1_0(29)) + (occ_func_1_0(46)*occ_func_1_0(5)) + (occ_func_1_0(30)*occ_func_1_0(2)) + (occ_func_1_0(50)*occ_func_1_0(29)) + (occ_func_1_0(39)*occ_func_1_0(1)) + (occ_func_1_0(30)*occ_func_1_0(6)) + (occ_func_1_0(33)*occ_func_1_0(27)) + (occ_func_1_0(56)*occ_func_1_0(5)) + (occ_func_1_0(32)*occ_func_1_0(2)) + (occ_func_1_0(47)*occ_func_1_0(31)) + (occ_func_1_0(42)*occ_func_1_0(4)) + (occ_func_1_0(28)*occ_func_1_0(3)) + (occ_func_1_0(38)*occ_func_1_0(29)) + (occ_func_1_0(51)*occ_func_1_0(6)) + (occ_func_1_0(30)*occ_func_1_0(1)) + (occ_func_1_0(53)*occ_func_1_0(32)) + (occ_func_1_0(36)*occ_func_1_0(5)) + (occ_func_1_0(27)*occ_func_1_0(2)) + (occ_func_1_0(42)*occ_func_1_0(28)) + (occ_func_1_0(47)*occ_func_1_0(3)) + (occ_func_1_0(31)*occ_func_1_0(4)) + (occ_func_1_0(54)*occ_func_1_0(32)) + (occ_func_1_0(35)*occ_func_1_0(4)) + (occ_func_1_0(27)*occ_func_1_0(3)) + (occ_func_1_0(52)*occ_func_1_0(31)) + (occ_func_1_0(37)*occ_func_1_0(1)) + (occ_func_1_0(28)*occ_func_1_0(6)) + (occ_func_1_0(40)*occ_func_1_0(31)) + (occ_func_1_0(49)*occ_func_1_0(6)) + (occ_func_1_0(28)*occ_func_1_0(1)) + (occ_func_1_0(34)*occ_func_1_0(27)) + (occ_func_1_0(55)*occ_func_1_0(4)) + (occ_func_1_0(32)*occ_func_1_0(3)) + (occ_func_1_0(49)*occ_func_1_0(28)) + (occ_func_1_0(40)*occ_func_1_0(1)) + (occ_func_1_0(31)*occ_func_1_0(6)) + (occ_func_1_0(35)*occ_func_1_0(27)) + (occ_func_1_0(54)*occ_func_1_0(3)) + (occ_func_1_0(32)*occ_func_1_0(4)) + (occ_func_1_0(55)*occ_func_1_0(32)) + (occ_func_1_0(34)*occ_func_1_0(3)) + (occ_func_1_0(27)*occ_func_1_0(4)) + (occ_func_1_0(37)*occ_func_1_0(28)) + (occ_func_1_0(52)*occ_func_1_0(6)) + (occ_func_1_0(31)*occ_func_1_0(1)) + (occ_func_1_0(45)*occ_func_1_0(29)) + (occ_func_1_0(44)*occ_func_1_0(2)) + (occ_func_1_0(30)*occ_func_1_0(5)) + (occ_func_1_0(44)*occ_func_1_0(30)) + (occ_func_1_0(45)*occ_func_1_0(5)) + (occ_func_1_0(29)*occ_func_1_0(2)) + (occ_func_1_0(46)*occ_func_1_0(30)) + (occ_func_1_0(43)*occ_func_1_0(2)) + (occ_func_1_0(29)*occ_func_1_0(5)) + (occ_func_1_0(48)*occ_func_1_0(31)) + (occ_func_1_0(41)*occ_func_1_0(3)) + (occ_func_1_0(28)*occ_func_1_0(4)) + (occ_func_1_0(41)*occ_func_1_0(28)) + (occ_func_1_0(48)*occ_func_1_0(4)) + (occ_func_1_0(31)*occ_func_1_0(3)) + (occ_func_1_0(51)*occ_func_1_0(30)) + (occ_func_1_0(38)*occ_func_1_0(1)) + (occ_func_1_0(29)*occ_func_1_0(6)) + (occ_func_1_0(56)*occ_func_1_0(32)) + (occ_func_1_0(33)*occ_func_1_0(2)) + (occ_func_1_0(27)*occ_func_1_0(5)) + (occ_func_1_0(36)*occ_func_1_0(27)) + (occ_func_1_0(53)*occ_func_1_0(2)) + (occ_func_1_0(32)*occ_func_1_0(5)) + (occ_func_1_0(39)*occ_func_1_0(30)) + (occ_func_1_0(50)*occ_func_1_0(6)) + (occ_func_1_0(29)*occ_func_1_0(1)))/24.0;
  }

  /**** Basis functions for orbit 3, 7****
#Points: 3
MaxLength: 8.8104314  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_7_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(15)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(18)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(15)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(9)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(17)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(14)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(14)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(17)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(18)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(9)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_7_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(11)) + (occ_func_1_0(46)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(14)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(16)) + (occ_func_1_0(39)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(9)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(7)) + (occ_func_1_0(56)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(18)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(13)) + (occ_func_1_0(42)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(12)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(8)) + (occ_func_1_0(51)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(17)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(15)) + (occ_func_1_0(36)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(10)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(12)) + (occ_func_1_0(47)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(13)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(16)) + (occ_func_1_0(35)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(9)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(18)) + (occ_func_1_0(37)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(7)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(10)) + (occ_func_1_0(49)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(15)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(8)) + (occ_func_1_0(55)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(17)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(15)) + (occ_func_1_0(40)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(10)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(9)) + (occ_func_1_0(54)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(16)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(17)) + (occ_func_1_0(34)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(8)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(7)) + (occ_func_1_0(52)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(18)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(13)) + (occ_func_1_0(44)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(12)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(12)) + (occ_func_1_0(45)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(13)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(14)) + (occ_func_1_0(43)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(11)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(14)) + (occ_func_1_0(41)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(11)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(11)) + (occ_func_1_0(48)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(14)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(17)) + (occ_func_1_0(38)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(8)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(18)) + (occ_func_1_0(33)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(7)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(10)) + (occ_func_1_0(53)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(15)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(9)) + (occ_func_1_0(50)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(16)*occ_func_1_0(4)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_7_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(43)*occ_func_1_0(11)) + (occ_func_1_0(46)*occ_func_1_0(4)) + (occ_func_1_0(14)*occ_func_1_0(3)) + (occ_func_1_0(50)*occ_func_1_0(16)) + (occ_func_1_0(39)*occ_func_1_0(4)) + (occ_func_1_0(9)*occ_func_1_0(3)) + (occ_func_1_0(33)*occ_func_1_0(7)) + (occ_func_1_0(56)*occ_func_1_0(6)) + (occ_func_1_0(18)*occ_func_1_0(1)) + (occ_func_1_0(47)*occ_func_1_0(13)) + (occ_func_1_0(42)*occ_func_1_0(2)) + (occ_func_1_0(12)*occ_func_1_0(5)) + (occ_func_1_0(38)*occ_func_1_0(8)) + (occ_func_1_0(51)*occ_func_1_0(4)) + (occ_func_1_0(17)*occ_func_1_0(3)) + (occ_func_1_0(53)*occ_func_1_0(15)) + (occ_func_1_0(36)*occ_func_1_0(1)) + (occ_func_1_0(10)*occ_func_1_0(6)) + (occ_func_1_0(42)*occ_func_1_0(12)) + (occ_func_1_0(47)*occ_func_1_0(5)) + (occ_func_1_0(13)*occ_func_1_0(2)) + (occ_func_1_0(54)*occ_func_1_0(16)) + (occ_func_1_0(35)*occ_func_1_0(1)) + (occ_func_1_0(9)*occ_func_1_0(6)) + (occ_func_1_0(52)*occ_func_1_0(18)) + (occ_func_1_0(37)*occ_func_1_0(2)) + (occ_func_1_0(7)*occ_func_1_0(5)) + (occ_func_1_0(40)*occ_func_1_0(10)) + (occ_func_1_0(49)*occ_func_1_0(2)) + (occ_func_1_0(15)*occ_func_1_0(5)) + (occ_func_1_0(34)*occ_func_1_0(8)) + (occ_func_1_0(55)*occ_func_1_0(6)) + (occ_func_1_0(17)*occ_func_1_0(1)) + (occ_func_1_0(49)*occ_func_1_0(15)) + (occ_func_1_0(40)*occ_func_1_0(5)) + (occ_func_1_0(10)*occ_func_1_0(2)) + (occ_func_1_0(35)*occ_func_1_0(9)) + (occ_func_1_0(54)*occ_func_1_0(6)) + (occ_func_1_0(16)*occ_func_1_0(1)) + (occ_func_1_0(55)*occ_func_1_0(17)) + (occ_func_1_0(34)*occ_func_1_0(1)) + (occ_func_1_0(8)*occ_func_1_0(6)) + (occ_func_1_0(37)*occ_func_1_0(7)) + (occ_func_1_0(52)*occ_func_1_0(5)) + (occ_func_1_0(18)*occ_func_1_0(2)) + (occ_func_1_0(45)*occ_func_1_0(13)) + (occ_func_1_0(44)*occ_func_1_0(4)) + (occ_func_1_0(12)*occ_func_1_0(3)) + (occ_func_1_0(44)*occ_func_1_0(12)) + (occ_func_1_0(45)*occ_func_1_0(3)) + (occ_func_1_0(13)*occ_func_1_0(4)) + (occ_func_1_0(46)*occ_func_1_0(14)) + (occ_func_1_0(43)*occ_func_1_0(3)) + (occ_func_1_0(11)*occ_func_1_0(4)) + (occ_func_1_0(48)*occ_func_1_0(14)) + (occ_func_1_0(41)*occ_func_1_0(2)) + (occ_func_1_0(11)*occ_func_1_0(5)) + (occ_func_1_0(41)*occ_func_1_0(11)) + (occ_func_1_0(48)*occ_func_1_0(5)) + (occ_func_1_0(14)*occ_func_1_0(2)) + (occ_func_1_0(51)*occ_func_1_0(17)) + (occ_func_1_0(38)*occ_func_1_0(3)) + (occ_func_1_0(8)*occ_func_1_0(4)) + (occ_func_1_0(56)*occ_func_1_0(18)) + (occ_func_1_0(33)*occ_func_1_0(1)) + (occ_func_1_0(7)*occ_func_1_0(6)) + (occ_func_1_0(36)*occ_func_1_0(10)) + (occ_func_1_0(53)*occ_func_1_0(6)) + (occ_func_1_0(15)*occ_func_1_0(1)) + (occ_func_1_0(39)*occ_func_1_0(9)) + (occ_func_1_0(50)*occ_func_1_0(3)) + (occ_func_1_0(16)*occ_func_1_0(4)))/24.0;
  }

  /**** Basis functions for orbit 3, 8****
#Points: 3
MaxLength: 8.8104314  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000   -1.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_8_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(44)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(39)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(48)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_8_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(41)) + (occ_func_1_0(46)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(48)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(54)) + (occ_func_1_0(39)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(35)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(37)) + (occ_func_1_0(56)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(52)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(45)) + (occ_func_1_0(42)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(44)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(34)) + (occ_func_1_0(51)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(55)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(49)) + (occ_func_1_0(36)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(40)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(44)) + (occ_func_1_0(47)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(45)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(56)) + (occ_func_1_0(37)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(33)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(36)) + (occ_func_1_0(49)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(53)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(39)) + (occ_func_1_0(54)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(50)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(51)) + (occ_func_1_0(34)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(38)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(48)) + (occ_func_1_0(43)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(41)*occ_func_1_0(12)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_8_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(43)*occ_func_1_0(41)) + (occ_func_1_0(46)*occ_func_1_0(12)) + (occ_func_1_0(48)*occ_func_1_0(13)) + (occ_func_1_0(50)*occ_func_1_0(54)) + (occ_func_1_0(39)*occ_func_1_0(17)) + (occ_func_1_0(35)*occ_func_1_0(8)) + (occ_func_1_0(33)*occ_func_1_0(37)) + (occ_func_1_0(56)*occ_func_1_0(15)) + (occ_func_1_0(52)*occ_func_1_0(10)) + (occ_func_1_0(47)*occ_func_1_0(45)) + (occ_func_1_0(42)*occ_func_1_0(11)) + (occ_func_1_0(44)*occ_func_1_0(14)) + (occ_func_1_0(38)*occ_func_1_0(34)) + (occ_func_1_0(51)*occ_func_1_0(9)) + (occ_func_1_0(55)*occ_func_1_0(16)) + (occ_func_1_0(53)*occ_func_1_0(49)) + (occ_func_1_0(36)*occ_func_1_0(7)) + (occ_func_1_0(40)*occ_func_1_0(18)) + (occ_func_1_0(42)*occ_func_1_0(44)) + (occ_func_1_0(47)*occ_func_1_0(14)) + (occ_func_1_0(45)*occ_func_1_0(11)) + (occ_func_1_0(52)*occ_func_1_0(56)) + (occ_func_1_0(37)*occ_func_1_0(15)) + (occ_func_1_0(33)*occ_func_1_0(10)) + (occ_func_1_0(40)*occ_func_1_0(36)) + (occ_func_1_0(49)*occ_func_1_0(7)) + (occ_func_1_0(53)*occ_func_1_0(18)) + (occ_func_1_0(35)*occ_func_1_0(39)) + (occ_func_1_0(54)*occ_func_1_0(17)) + (occ_func_1_0(50)*occ_func_1_0(8)) + (occ_func_1_0(55)*occ_func_1_0(51)) + (occ_func_1_0(34)*occ_func_1_0(9)) + (occ_func_1_0(38)*occ_func_1_0(16)) + (occ_func_1_0(46)*occ_func_1_0(48)) + (occ_func_1_0(43)*occ_func_1_0(13)) + (occ_func_1_0(41)*occ_func_1_0(12)))/12.0;
  }

  /**** Basis functions for orbit 3, 9****
#Points: 3
MaxLength: 8.8104314  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
              -0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_9_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(26)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_9_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(19)) + (occ_func_1_0(46)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(26)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(23)) + (occ_func_1_0(39)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(22)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(20)) + (occ_func_1_0(56)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(25)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(21)) + (occ_func_1_0(42)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(24)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(21)) + (occ_func_1_0(51)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(24)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(23)) + (occ_func_1_0(36)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(22)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(20)) + (occ_func_1_0(47)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(25)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(25)) + (occ_func_1_0(35)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(20)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(25)) + (occ_func_1_0(37)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(20)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(22)) + (occ_func_1_0(49)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(23)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(19)) + (occ_func_1_0(55)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(26)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(24)) + (occ_func_1_0(40)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(21)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(22)) + (occ_func_1_0(54)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(23)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(24)) + (occ_func_1_0(34)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(21)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(19)) + (occ_func_1_0(52)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(26)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(25)) + (occ_func_1_0(44)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(20)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(24)) + (occ_func_1_0(45)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(21)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(22)) + (occ_func_1_0(43)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(23)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(26)) + (occ_func_1_0(41)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(19)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(23)) + (occ_func_1_0(48)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(22)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(26)) + (occ_func_1_0(38)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(19)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(26)) + (occ_func_1_0(33)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(19)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(21)) + (occ_func_1_0(53)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(24)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(20)) + (occ_func_1_0(50)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(25)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(20)) + (occ_func_1_0(45)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(25)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(21)) + (occ_func_1_0(44)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(24)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(23)) + (occ_func_1_0(46)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(22)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(19)) + (occ_func_1_0(48)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(26)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(22)) + (occ_func_1_0(41)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(23)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(19)) + (occ_func_1_0(51)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(26)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(19)) + (occ_func_1_0(56)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(26)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(24)) + (occ_func_1_0(36)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(21)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(25)) + (occ_func_1_0(39)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(20)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(20)) + (occ_func_1_0(54)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(25)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(20)) + (occ_func_1_0(52)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(25)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(23)) + (occ_func_1_0(40)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(22)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(26)) + (occ_func_1_0(34)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(19)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(21)) + (occ_func_1_0(49)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(24)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(23)) + (occ_func_1_0(35)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(22)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(21)) + (occ_func_1_0(55)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(24)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(26)) + (occ_func_1_0(37)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(19)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(22)) + (occ_func_1_0(50)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(23)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(25)) + (occ_func_1_0(33)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(20)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(24)) + (occ_func_1_0(47)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(21)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(24)) + (occ_func_1_0(38)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(21)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(22)) + (occ_func_1_0(53)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(23)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(25)) + (occ_func_1_0(42)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(20)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(26)) + (occ_func_1_0(43)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(19)*occ_func_1_0(9)*occ_func_1_0(0)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_9_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(43)*occ_func_1_0(19)) + (occ_func_1_0(46)*occ_func_1_0(9)) + (occ_func_1_0(26)*occ_func_1_0(16)) + (occ_func_1_0(50)*occ_func_1_0(23)) + (occ_func_1_0(39)*occ_func_1_0(12)) + (occ_func_1_0(22)*occ_func_1_0(13)) + (occ_func_1_0(33)*occ_func_1_0(20)) + (occ_func_1_0(56)*occ_func_1_0(17)) + (occ_func_1_0(25)*occ_func_1_0(8)) + (occ_func_1_0(47)*occ_func_1_0(21)) + (occ_func_1_0(42)*occ_func_1_0(7)) + (occ_func_1_0(24)*occ_func_1_0(18)) + (occ_func_1_0(38)*occ_func_1_0(21)) + (occ_func_1_0(51)*occ_func_1_0(14)) + (occ_func_1_0(24)*occ_func_1_0(11)) + (occ_func_1_0(53)*occ_func_1_0(23)) + (occ_func_1_0(36)*occ_func_1_0(8)) + (occ_func_1_0(22)*occ_func_1_0(17)) + (occ_func_1_0(42)*occ_func_1_0(20)) + (occ_func_1_0(47)*occ_func_1_0(10)) + (occ_func_1_0(25)*occ_func_1_0(15)) + (occ_func_1_0(54)*occ_func_1_0(25)) + (occ_func_1_0(35)*occ_func_1_0(10)) + (occ_func_1_0(20)*occ_func_1_0(15)) + (occ_func_1_0(52)*occ_func_1_0(25)) + (occ_func_1_0(37)*occ_func_1_0(11)) + (occ_func_1_0(20)*occ_func_1_0(14)) + (occ_func_1_0(40)*occ_func_1_0(22)) + (occ_func_1_0(49)*occ_func_1_0(12)) + (occ_func_1_0(23)*occ_func_1_0(13)) + (occ_func_1_0(34)*occ_func_1_0(19)) + (occ_func_1_0(55)*occ_func_1_0(15)) + (occ_func_1_0(26)*occ_func_1_0(10)) + (occ_func_1_0(49)*occ_func_1_0(24)) + (occ_func_1_0(40)*occ_func_1_0(14)) + (occ_func_1_0(21)*occ_func_1_0(11)) + (occ_func_1_0(35)*occ_func_1_0(22)) + (occ_func_1_0(54)*occ_func_1_0(18)) + (occ_func_1_0(23)*occ_func_1_0(7)) + (occ_func_1_0(55)*occ_func_1_0(24)) + (occ_func_1_0(34)*occ_func_1_0(7)) + (occ_func_1_0(21)*occ_func_1_0(18)) + (occ_func_1_0(37)*occ_func_1_0(19)) + (occ_func_1_0(52)*occ_func_1_0(13)) + (occ_func_1_0(26)*occ_func_1_0(12)) + (occ_func_1_0(45)*occ_func_1_0(25)) + (occ_func_1_0(44)*occ_func_1_0(17)) + (occ_func_1_0(20)*occ_func_1_0(8)) + (occ_func_1_0(44)*occ_func_1_0(24)) + (occ_func_1_0(45)*occ_func_1_0(16)) + (occ_func_1_0(21)*occ_func_1_0(9)) + (occ_func_1_0(46)*occ_func_1_0(22)) + (occ_func_1_0(43)*occ_func_1_0(8)) + (occ_func_1_0(23)*occ_func_1_0(17)) + (occ_func_1_0(48)*occ_func_1_0(26)) + (occ_func_1_0(41)*occ_func_1_0(15)) + (occ_func_1_0(19)*occ_func_1_0(10)) + (occ_func_1_0(41)*occ_func_1_0(23)) + (occ_func_1_0(48)*occ_func_1_0(18)) + (occ_func_1_0(22)*occ_func_1_0(7)) + (occ_func_1_0(51)*occ_func_1_0(26)) + (occ_func_1_0(38)*occ_func_1_0(13)) + (occ_func_1_0(19)*occ_func_1_0(12)) + (occ_func_1_0(56)*occ_func_1_0(26)) + (occ_func_1_0(33)*occ_func_1_0(9)) + (occ_func_1_0(19)*occ_func_1_0(16)) + (occ_func_1_0(36)*occ_func_1_0(21)) + (occ_func_1_0(53)*occ_func_1_0(16)) + (occ_func_1_0(24)*occ_func_1_0(9)) + (occ_func_1_0(39)*occ_func_1_0(20)) + (occ_func_1_0(50)*occ_func_1_0(11)) + (occ_func_1_0(25)*occ_func_1_0(14)) + (occ_func_1_0(44)*occ_func_1_0(20)) + (occ_func_1_0(45)*occ_func_1_0(8)) + (occ_func_1_0(25)*occ_func_1_0(17)) + (occ_func_1_0(45)*occ_func_1_0(21)) + (occ_func_1_0(44)*occ_func_1_0(9)) + (occ_func_1_0(24)*occ_func_1_0(16)) + (occ_func_1_0(43)*occ_func_1_0(23)) + (occ_func_1_0(46)*occ_func_1_0(17)) + (occ_func_1_0(22)*occ_func_1_0(8)) + (occ_func_1_0(41)*occ_func_1_0(19)) + (occ_func_1_0(48)*occ_func_1_0(10)) + (occ_func_1_0(26)*occ_func_1_0(15)) + (occ_func_1_0(48)*occ_func_1_0(22)) + (occ_func_1_0(41)*occ_func_1_0(7)) + (occ_func_1_0(23)*occ_func_1_0(18)) + (occ_func_1_0(38)*occ_func_1_0(19)) + (occ_func_1_0(51)*occ_func_1_0(12)) + (occ_func_1_0(26)*occ_func_1_0(13)) + (occ_func_1_0(33)*occ_func_1_0(19)) + (occ_func_1_0(56)*occ_func_1_0(16)) + (occ_func_1_0(26)*occ_func_1_0(9)) + (occ_func_1_0(53)*occ_func_1_0(24)) + (occ_func_1_0(36)*occ_func_1_0(9)) + (occ_func_1_0(21)*occ_func_1_0(16)) + (occ_func_1_0(50)*occ_func_1_0(25)) + (occ_func_1_0(39)*occ_func_1_0(14)) + (occ_func_1_0(20)*occ_func_1_0(11)) + (occ_func_1_0(35)*occ_func_1_0(20)) + (occ_func_1_0(54)*occ_func_1_0(15)) + (occ_func_1_0(25)*occ_func_1_0(10)) + (occ_func_1_0(37)*occ_func_1_0(20)) + (occ_func_1_0(52)*occ_func_1_0(14)) + (occ_func_1_0(25)*occ_func_1_0(11)) + (occ_func_1_0(49)*occ_func_1_0(23)) + (occ_func_1_0(40)*occ_func_1_0(13)) + (occ_func_1_0(22)*occ_func_1_0(12)) + (occ_func_1_0(55)*occ_func_1_0(26)) + (occ_func_1_0(34)*occ_func_1_0(10)) + (occ_func_1_0(19)*occ_func_1_0(15)) + (occ_func_1_0(40)*occ_func_1_0(21)) + (occ_func_1_0(49)*occ_func_1_0(11)) + (occ_func_1_0(24)*occ_func_1_0(14)) + (occ_func_1_0(54)*occ_func_1_0(23)) + (occ_func_1_0(35)*occ_func_1_0(7)) + (occ_func_1_0(22)*occ_func_1_0(18)) + (occ_func_1_0(34)*occ_func_1_0(21)) + (occ_func_1_0(55)*occ_func_1_0(18)) + (occ_func_1_0(24)*occ_func_1_0(7)) + (occ_func_1_0(52)*occ_func_1_0(26)) + (occ_func_1_0(37)*occ_func_1_0(12)) + (occ_func_1_0(19)*occ_func_1_0(13)) + (occ_func_1_0(39)*occ_func_1_0(22)) + (occ_func_1_0(50)*occ_func_1_0(13)) + (occ_func_1_0(23)*occ_func_1_0(12)) + (occ_func_1_0(56)*occ_func_1_0(25)) + (occ_func_1_0(33)*occ_func_1_0(8)) + (occ_func_1_0(20)*occ_func_1_0(17)) + (occ_func_1_0(42)*occ_func_1_0(24)) + (occ_func_1_0(47)*occ_func_1_0(18)) + (occ_func_1_0(21)*occ_func_1_0(7)) + (occ_func_1_0(51)*occ_func_1_0(24)) + (occ_func_1_0(38)*occ_func_1_0(11)) + (occ_func_1_0(21)*occ_func_1_0(14)) + (occ_func_1_0(36)*occ_func_1_0(22)) + (occ_func_1_0(53)*occ_func_1_0(17)) + (occ_func_1_0(23)*occ_func_1_0(8)) + (occ_func_1_0(47)*occ_func_1_0(25)) + (occ_func_1_0(42)*occ_func_1_0(15)) + (occ_func_1_0(20)*occ_func_1_0(10)) + (occ_func_1_0(46)*occ_func_1_0(26)) + (occ_func_1_0(43)*occ_func_1_0(16)) + (occ_func_1_0(19)*occ_func_1_0(9)))/48.0;
  }

  /**** Basis functions for orbit 3, 10****
#Points: 3
MaxLength: 8.8104314  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
              -0.5000000    0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_10_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(48)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(42)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(39)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(46)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(44)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_10_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(38)) + (occ_func_1_0(46)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(51)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(43)) + (occ_func_1_0(39)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(46)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(35)) + (occ_func_1_0(56)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(54)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(40)) + (occ_func_1_0(42)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(49)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(38)*occ_func_1_0(45)) + (occ_func_1_0(51)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(44)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(54)) + (occ_func_1_0(36)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(35)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(37)) + (occ_func_1_0(47)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(52)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(56)) + (occ_func_1_0(35)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(33)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(47)) + (occ_func_1_0(37)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(42)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(40)*occ_func_1_0(48)) + (occ_func_1_0(49)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(41)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(33)) + (occ_func_1_0(55)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(56)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(42)) + (occ_func_1_0(40)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(47)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(35)*occ_func_1_0(36)) + (occ_func_1_0(54)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(53)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(55)*occ_func_1_0(53)) + (occ_func_1_0(34)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(36)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(37)*occ_func_1_0(41)) + (occ_func_1_0(52)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(48)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(45)*occ_func_1_0(50)) + (occ_func_1_0(44)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(39)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(51)) + (occ_func_1_0(45)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(38)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(46)*occ_func_1_0(39)) + (occ_func_1_0(43)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(50)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(48)*occ_func_1_0(52)) + (occ_func_1_0(41)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(37)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(41)*occ_func_1_0(49)) + (occ_func_1_0(48)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(40)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(46)) + (occ_func_1_0(38)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(43)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(56)*occ_func_1_0(55)) + (occ_func_1_0(33)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(34)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(36)*occ_func_1_0(34)) + (occ_func_1_0(53)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(55)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(39)*occ_func_1_0(44)) + (occ_func_1_0(50)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(45)*occ_func_1_0(10)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_10_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(43)*occ_func_1_0(38)) + (occ_func_1_0(46)*occ_func_1_0(10)) + (occ_func_1_0(51)*occ_func_1_0(15)) + (occ_func_1_0(50)*occ_func_1_0(43)) + (occ_func_1_0(39)*occ_func_1_0(7)) + (occ_func_1_0(46)*occ_func_1_0(18)) + (occ_func_1_0(33)*occ_func_1_0(35)) + (occ_func_1_0(56)*occ_func_1_0(14)) + (occ_func_1_0(54)*occ_func_1_0(11)) + (occ_func_1_0(47)*occ_func_1_0(40)) + (occ_func_1_0(42)*occ_func_1_0(9)) + (occ_func_1_0(49)*occ_func_1_0(16)) + (occ_func_1_0(38)*occ_func_1_0(45)) + (occ_func_1_0(51)*occ_func_1_0(18)) + (occ_func_1_0(44)*occ_func_1_0(7)) + (occ_func_1_0(53)*occ_func_1_0(54)) + (occ_func_1_0(36)*occ_func_1_0(13)) + (occ_func_1_0(35)*occ_func_1_0(12)) + (occ_func_1_0(42)*occ_func_1_0(37)) + (occ_func_1_0(47)*occ_func_1_0(8)) + (occ_func_1_0(52)*occ_func_1_0(17)) + (occ_func_1_0(54)*occ_func_1_0(56)) + (occ_func_1_0(35)*occ_func_1_0(14)) + (occ_func_1_0(33)*occ_func_1_0(11)) + (occ_func_1_0(52)*occ_func_1_0(47)) + (occ_func_1_0(37)*occ_func_1_0(8)) + (occ_func_1_0(42)*occ_func_1_0(17)) + (occ_func_1_0(40)*occ_func_1_0(48)) + (occ_func_1_0(49)*occ_func_1_0(17)) + (occ_func_1_0(41)*occ_func_1_0(8)) + (occ_func_1_0(34)*occ_func_1_0(33)) + (occ_func_1_0(55)*occ_func_1_0(12)) + (occ_func_1_0(56)*occ_func_1_0(13)) + (occ_func_1_0(49)*occ_func_1_0(42)) + (occ_func_1_0(40)*occ_func_1_0(9)) + (occ_func_1_0(47)*occ_func_1_0(16)) + (occ_func_1_0(35)*occ_func_1_0(36)) + (occ_func_1_0(54)*occ_func_1_0(13)) + (occ_func_1_0(53)*occ_func_1_0(12)) + (occ_func_1_0(55)*occ_func_1_0(53)) + (occ_func_1_0(34)*occ_func_1_0(11)) + (occ_func_1_0(36)*occ_func_1_0(14)) + (occ_func_1_0(37)*occ_func_1_0(41)) + (occ_func_1_0(52)*occ_func_1_0(16)) + (occ_func_1_0(48)*occ_func_1_0(9)) + (occ_func_1_0(45)*occ_func_1_0(50)) + (occ_func_1_0(44)*occ_func_1_0(15)) + (occ_func_1_0(39)*occ_func_1_0(10)) + (occ_func_1_0(44)*occ_func_1_0(51)) + (occ_func_1_0(45)*occ_func_1_0(18)) + (occ_func_1_0(38)*occ_func_1_0(7)) + (occ_func_1_0(46)*occ_func_1_0(39)) + (occ_func_1_0(43)*occ_func_1_0(7)) + (occ_func_1_0(50)*occ_func_1_0(18)) + (occ_func_1_0(48)*occ_func_1_0(52)) + (occ_func_1_0(41)*occ_func_1_0(16)) + (occ_func_1_0(37)*occ_func_1_0(9)) + (occ_func_1_0(41)*occ_func_1_0(49)) + (occ_func_1_0(48)*occ_func_1_0(17)) + (occ_func_1_0(40)*occ_func_1_0(8)) + (occ_func_1_0(51)*occ_func_1_0(46)) + (occ_func_1_0(38)*occ_func_1_0(10)) + (occ_func_1_0(43)*occ_func_1_0(15)) + (occ_func_1_0(56)*occ_func_1_0(55)) + (occ_func_1_0(33)*occ_func_1_0(12)) + (occ_func_1_0(34)*occ_func_1_0(13)) + (occ_func_1_0(36)*occ_func_1_0(34)) + (occ_func_1_0(53)*occ_func_1_0(11)) + (occ_func_1_0(55)*occ_func_1_0(14)) + (occ_func_1_0(39)*occ_func_1_0(44)) + (occ_func_1_0(50)*occ_func_1_0(15)) + (occ_func_1_0(45)*occ_func_1_0(10)))/24.0;
  }

  /**** Basis functions for orbit 3, 11****
#Points: 3
MaxLength: 8.8104314  MinLength: 7.8802894
               0.5000000    0.5000000    0.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000    1.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_11_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(48)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(46)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(39)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_11_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(43)*occ_func_1_0(45)) + (occ_func_1_0(46)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(44)*occ_func_1_0(28)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(50)*occ_func_1_0(38)) + (occ_func_1_0(39)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(51)*occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(33)*occ_func_1_0(36)) + (occ_func_1_0(56)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(53)*occ_func_1_0(28)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(47)*occ_func_1_0(48)) + (occ_func_1_0(42)*occ_func_1_0(0)*occ_func_1_0(30)) + (occ_func_1_0(41)*occ_func_1_0(29)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(53)*occ_func_1_0(56)) + (occ_func_1_0(36)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(33)*occ_func_1_0(28)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(42)*occ_func_1_0(41)) + (occ_func_1_0(47)*occ_func_1_0(0)*occ_func_1_0(29)) + (occ_func_1_0(48)*occ_func_1_0(30)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(54)*occ_func_1_0(55)) + (occ_func_1_0(35)*occ_func_1_0(0)*occ_func_1_0(30)) + (occ_func_1_0(34)*occ_func_1_0(29)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(52)*occ_func_1_0(40)) + (occ_func_1_0(37)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(49)*occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(34)*occ_func_1_0(35)) + (occ_func_1_0(55)*occ_func_1_0(0)*occ_func_1_0(30)) + (occ_func_1_0(54)*occ_func_1_0(29)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(49)*occ_func_1_0(37)) + (occ_func_1_0(40)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(52)*occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(44)*occ_func_1_0(46)) + (occ_func_1_0(45)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(43)*occ_func_1_0(28)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(51)*occ_func_1_0(39)) + (occ_func_1_0(38)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(50)*occ_func_1_0(32)*occ_func_1_0(0)))/12.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_11_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(43)*occ_func_1_0(45)) + (occ_func_1_0(46)*occ_func_1_0(31)) + (occ_func_1_0(44)*occ_func_1_0(28)) + (occ_func_1_0(50)*occ_func_1_0(38)) + (occ_func_1_0(39)*occ_func_1_0(27)) + (occ_func_1_0(51)*occ_func_1_0(32)) + (occ_func_1_0(33)*occ_func_1_0(36)) + (occ_func_1_0(56)*occ_func_1_0(31)) + (occ_func_1_0(53)*occ_func_1_0(28)) + (occ_func_1_0(47)*occ_func_1_0(48)) + (occ_func_1_0(42)*occ_func_1_0(30)) + (occ_func_1_0(41)*occ_func_1_0(29)) + (occ_func_1_0(53)*occ_func_1_0(56)) + (occ_func_1_0(36)*occ_func_1_0(31)) + (occ_func_1_0(33)*occ_func_1_0(28)) + (occ_func_1_0(42)*occ_func_1_0(41)) + (occ_func_1_0(47)*occ_func_1_0(29)) + (occ_func_1_0(48)*occ_func_1_0(30)) + (occ_func_1_0(54)*occ_func_1_0(55)) + (occ_func_1_0(35)*occ_func_1_0(30)) + (occ_func_1_0(34)*occ_func_1_0(29)) + (occ_func_1_0(52)*occ_func_1_0(40)) + (occ_func_1_0(37)*occ_func_1_0(27)) + (occ_func_1_0(49)*occ_func_1_0(32)) + (occ_func_1_0(34)*occ_func_1_0(35)) + (occ_func_1_0(55)*occ_func_1_0(30)) + (occ_func_1_0(54)*occ_func_1_0(29)) + (occ_func_1_0(49)*occ_func_1_0(37)) + (occ_func_1_0(40)*occ_func_1_0(27)) + (occ_func_1_0(52)*occ_func_1_0(32)) + (occ_func_1_0(44)*occ_func_1_0(46)) + (occ_func_1_0(45)*occ_func_1_0(31)) + (occ_func_1_0(43)*occ_func_1_0(28)) + (occ_func_1_0(51)*occ_func_1_0(39)) + (occ_func_1_0(38)*occ_func_1_0(27)) + (occ_func_1_0(50)*occ_func_1_0(32)))/12.0;
  }

  /**** Basis functions for orbit 3, 12****
#Points: 3
MaxLength: 9.6513440  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000   -0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_12_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(42)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(44)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(46)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(48)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(39)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(44)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(48)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(39)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(42)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(46)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_12_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(43)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(46)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(50)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(39)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(33)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(56)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(47)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(42)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(38)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(51)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(53)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(36)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(42)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(47)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(54)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(35)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(52)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(37)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(40)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(49)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(34)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(55)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(49)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(40)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(35)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(54)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(55)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(34)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(37)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(52)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(45)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(44)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(44)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(45)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(46)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(43)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(48)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(41)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(41)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(48)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(51)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(38)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(56)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(33)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(36)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(53)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(39)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(50)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(44)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(45)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(45)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(44)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(43)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(46)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(41)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(48)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(48)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(41)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(38)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(51)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(33)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(56)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(53)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(36)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(50)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(39)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(35)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(54)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(37)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(52)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(49)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(40)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(55)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(34)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(40)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(49)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(54)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(35)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(34)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(55)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(52)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(37)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(39)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(50)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(56)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(33)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(42)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(47)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(51)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(38)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(36)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(53)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(47)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(42)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(46)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(43)*occ_func_1_0(6)*occ_func_1_0(0)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_12_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(43)) + (occ_func_1_0(74)*occ_func_1_0(6)) + (occ_func_1_0(46)*occ_func_1_0(1)) + (occ_func_1_0(71)*occ_func_1_0(50)) + (occ_func_1_0(66)*occ_func_1_0(5)) + (occ_func_1_0(39)*occ_func_1_0(2)) + (occ_func_1_0(58)*occ_func_1_0(33)) + (occ_func_1_0(79)*occ_func_1_0(3)) + (occ_func_1_0(56)*occ_func_1_0(4)) + (occ_func_1_0(67)*occ_func_1_0(47)) + (occ_func_1_0(70)*occ_func_1_0(6)) + (occ_func_1_0(42)*occ_func_1_0(1)) + (occ_func_1_0(65)*occ_func_1_0(38)) + (occ_func_1_0(72)*occ_func_1_0(2)) + (occ_func_1_0(51)*occ_func_1_0(5)) + (occ_func_1_0(77)*occ_func_1_0(53)) + (occ_func_1_0(60)*occ_func_1_0(4)) + (occ_func_1_0(36)*occ_func_1_0(3)) + (occ_func_1_0(62)*occ_func_1_0(42)) + (occ_func_1_0(75)*occ_func_1_0(6)) + (occ_func_1_0(47)*occ_func_1_0(1)) + (occ_func_1_0(79)*occ_func_1_0(54)) + (occ_func_1_0(58)*occ_func_1_0(2)) + (occ_func_1_0(35)*occ_func_1_0(5)) + (occ_func_1_0(75)*occ_func_1_0(52)) + (occ_func_1_0(62)*occ_func_1_0(4)) + (occ_func_1_0(37)*occ_func_1_0(3)) + (occ_func_1_0(68)*occ_func_1_0(40)) + (occ_func_1_0(69)*occ_func_1_0(3)) + (occ_func_1_0(49)*occ_func_1_0(4)) + (occ_func_1_0(57)*occ_func_1_0(34)) + (occ_func_1_0(80)*occ_func_1_0(5)) + (occ_func_1_0(55)*occ_func_1_0(2)) + (occ_func_1_0(70)*occ_func_1_0(49)) + (occ_func_1_0(67)*occ_func_1_0(3)) + (occ_func_1_0(40)*occ_func_1_0(4)) + (occ_func_1_0(60)*occ_func_1_0(35)) + (occ_func_1_0(77)*occ_func_1_0(2)) + (occ_func_1_0(54)*occ_func_1_0(5)) + (occ_func_1_0(78)*occ_func_1_0(55)) + (occ_func_1_0(59)*occ_func_1_0(5)) + (occ_func_1_0(34)*occ_func_1_0(2)) + (occ_func_1_0(61)*occ_func_1_0(37)) + (occ_func_1_0(76)*occ_func_1_0(4)) + (occ_func_1_0(52)*occ_func_1_0(3)) + (occ_func_1_0(73)*occ_func_1_0(45)) + (occ_func_1_0(64)*occ_func_1_0(1)) + (occ_func_1_0(44)*occ_func_1_0(6)) + (occ_func_1_0(72)*occ_func_1_0(44)) + (occ_func_1_0(65)*occ_func_1_0(1)) + (occ_func_1_0(45)*occ_func_1_0(6)) + (occ_func_1_0(66)*occ_func_1_0(46)) + (occ_func_1_0(71)*occ_func_1_0(6)) + (occ_func_1_0(43)*occ_func_1_0(1)) + (occ_func_1_0(76)*occ_func_1_0(48)) + (occ_func_1_0(61)*occ_func_1_0(1)) + (occ_func_1_0(41)*occ_func_1_0(6)) + (occ_func_1_0(69)*occ_func_1_0(41)) + (occ_func_1_0(68)*occ_func_1_0(1)) + (occ_func_1_0(48)*occ_func_1_0(6)) + (occ_func_1_0(74)*occ_func_1_0(51)) + (occ_func_1_0(63)*occ_func_1_0(2)) + (occ_func_1_0(38)*occ_func_1_0(5)) + (occ_func_1_0(80)*occ_func_1_0(56)) + (occ_func_1_0(57)*occ_func_1_0(3)) + (occ_func_1_0(33)*occ_func_1_0(4)) + (occ_func_1_0(59)*occ_func_1_0(36)) + (occ_func_1_0(78)*occ_func_1_0(4)) + (occ_func_1_0(53)*occ_func_1_0(3)) + (occ_func_1_0(64)*occ_func_1_0(39)) + (occ_func_1_0(73)*occ_func_1_0(5)) + (occ_func_1_0(50)*occ_func_1_0(2)) + (occ_func_1_0(64)*occ_func_1_0(44)) + (occ_func_1_0(73)*occ_func_1_0(6)) + (occ_func_1_0(45)*occ_func_1_0(1)) + (occ_func_1_0(65)*occ_func_1_0(45)) + (occ_func_1_0(72)*occ_func_1_0(6)) + (occ_func_1_0(44)*occ_func_1_0(1)) + (occ_func_1_0(71)*occ_func_1_0(43)) + (occ_func_1_0(66)*occ_func_1_0(1)) + (occ_func_1_0(46)*occ_func_1_0(6)) + (occ_func_1_0(61)*occ_func_1_0(41)) + (occ_func_1_0(76)*occ_func_1_0(6)) + (occ_func_1_0(48)*occ_func_1_0(1)) + (occ_func_1_0(68)*occ_func_1_0(48)) + (occ_func_1_0(69)*occ_func_1_0(6)) + (occ_func_1_0(41)*occ_func_1_0(1)) + (occ_func_1_0(63)*occ_func_1_0(38)) + (occ_func_1_0(74)*occ_func_1_0(5)) + (occ_func_1_0(51)*occ_func_1_0(2)) + (occ_func_1_0(57)*occ_func_1_0(33)) + (occ_func_1_0(80)*occ_func_1_0(4)) + (occ_func_1_0(56)*occ_func_1_0(3)) + (occ_func_1_0(78)*occ_func_1_0(53)) + (occ_func_1_0(59)*occ_func_1_0(3)) + (occ_func_1_0(36)*occ_func_1_0(4)) + (occ_func_1_0(73)*occ_func_1_0(50)) + (occ_func_1_0(64)*occ_func_1_0(2)) + (occ_func_1_0(39)*occ_func_1_0(5)) + (occ_func_1_0(58)*occ_func_1_0(35)) + (occ_func_1_0(79)*occ_func_1_0(5)) + (occ_func_1_0(54)*occ_func_1_0(2)) + (occ_func_1_0(62)*occ_func_1_0(37)) + (occ_func_1_0(75)*occ_func_1_0(3)) + (occ_func_1_0(52)*occ_func_1_0(4)) + (occ_func_1_0(69)*occ_func_1_0(49)) + (occ_func_1_0(68)*occ_func_1_0(4)) + (occ_func_1_0(40)*occ_func_1_0(3)) + (occ_func_1_0(80)*occ_func_1_0(55)) + (occ_func_1_0(57)*occ_func_1_0(2)) + (occ_func_1_0(34)*occ_func_1_0(5)) + (occ_func_1_0(67)*occ_func_1_0(40)) + (occ_func_1_0(70)*occ_func_1_0(4)) + (occ_func_1_0(49)*occ_func_1_0(3)) + (occ_func_1_0(77)*occ_func_1_0(54)) + (occ_func_1_0(60)*occ_func_1_0(5)) + (occ_func_1_0(35)*occ_func_1_0(2)) + (occ_func_1_0(59)*occ_func_1_0(34)) + (occ_func_1_0(78)*occ_func_1_0(2)) + (occ_func_1_0(55)*occ_func_1_0(5)) + (occ_func_1_0(76)*occ_func_1_0(52)) + (occ_func_1_0(61)*occ_func_1_0(3)) + (occ_func_1_0(37)*occ_func_1_0(4)) + (occ_func_1_0(66)*occ_func_1_0(39)) + (occ_func_1_0(71)*occ_func_1_0(2)) + (occ_func_1_0(50)*occ_func_1_0(5)) + (occ_func_1_0(79)*occ_func_1_0(56)) + (occ_func_1_0(58)*occ_func_1_0(4)) + (occ_func_1_0(33)*occ_func_1_0(3)) + (occ_func_1_0(70)*occ_func_1_0(42)) + (occ_func_1_0(67)*occ_func_1_0(1)) + (occ_func_1_0(47)*occ_func_1_0(6)) + (occ_func_1_0(72)*occ_func_1_0(51)) + (occ_func_1_0(65)*occ_func_1_0(5)) + (occ_func_1_0(38)*occ_func_1_0(2)) + (occ_func_1_0(60)*occ_func_1_0(36)) + (occ_func_1_0(77)*occ_func_1_0(3)) + (occ_func_1_0(53)*occ_func_1_0(4)) + (occ_func_1_0(75)*occ_func_1_0(47)) + (occ_func_1_0(62)*occ_func_1_0(1)) + (occ_func_1_0(42)*occ_func_1_0(6)) + (occ_func_1_0(74)*occ_func_1_0(46)) + (occ_func_1_0(63)*occ_func_1_0(1)) + (occ_func_1_0(43)*occ_func_1_0(6)))/48.0;
  }

  /**** Basis functions for orbit 3, 13****
#Points: 3
MaxLength: 9.6513440  MinLength: 3.9401447
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
              -0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_13_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(20)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(19)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(25)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(24)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(22)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(23)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(26)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(21)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(20)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_13_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(19)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(26)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(23)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(22)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(20)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(25)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(21)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(24)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(21)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(24)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(23)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(22)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(20)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(25)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(25)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(20)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(25)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(20)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(22)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(23)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(19)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(26)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(24)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(21)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(22)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(23)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(24)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(21)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(19)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(26)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(25)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(4)) + (occ_func_1_0(20)*occ_func_1_0(3)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(24)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(21)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(22)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(23)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(26)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(2)) + (occ_func_1_0(19)*occ_func_1_0(5)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(23)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(5)) + (occ_func_1_0(22)*occ_func_1_0(2)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(26)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(19)*occ_func_1_0(4)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(26)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(1)) + (occ_func_1_0(19)*occ_func_1_0(6)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(21)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(6)) + (occ_func_1_0(24)*occ_func_1_0(1)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(20)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(3)) + (occ_func_1_0(25)*occ_func_1_0(4)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_13_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(19)) + (occ_func_1_0(74)*occ_func_1_0(4)) + (occ_func_1_0(26)*occ_func_1_0(3)) + (occ_func_1_0(71)*occ_func_1_0(23)) + (occ_func_1_0(66)*occ_func_1_0(4)) + (occ_func_1_0(22)*occ_func_1_0(3)) + (occ_func_1_0(58)*occ_func_1_0(20)) + (occ_func_1_0(79)*occ_func_1_0(6)) + (occ_func_1_0(25)*occ_func_1_0(1)) + (occ_func_1_0(67)*occ_func_1_0(21)) + (occ_func_1_0(70)*occ_func_1_0(2)) + (occ_func_1_0(24)*occ_func_1_0(5)) + (occ_func_1_0(65)*occ_func_1_0(21)) + (occ_func_1_0(72)*occ_func_1_0(4)) + (occ_func_1_0(24)*occ_func_1_0(3)) + (occ_func_1_0(77)*occ_func_1_0(23)) + (occ_func_1_0(60)*occ_func_1_0(1)) + (occ_func_1_0(22)*occ_func_1_0(6)) + (occ_func_1_0(62)*occ_func_1_0(20)) + (occ_func_1_0(75)*occ_func_1_0(5)) + (occ_func_1_0(25)*occ_func_1_0(2)) + (occ_func_1_0(79)*occ_func_1_0(25)) + (occ_func_1_0(58)*occ_func_1_0(1)) + (occ_func_1_0(20)*occ_func_1_0(6)) + (occ_func_1_0(75)*occ_func_1_0(25)) + (occ_func_1_0(62)*occ_func_1_0(2)) + (occ_func_1_0(20)*occ_func_1_0(5)) + (occ_func_1_0(68)*occ_func_1_0(22)) + (occ_func_1_0(69)*occ_func_1_0(2)) + (occ_func_1_0(23)*occ_func_1_0(5)) + (occ_func_1_0(57)*occ_func_1_0(19)) + (occ_func_1_0(80)*occ_func_1_0(6)) + (occ_func_1_0(26)*occ_func_1_0(1)) + (occ_func_1_0(70)*occ_func_1_0(24)) + (occ_func_1_0(67)*occ_func_1_0(5)) + (occ_func_1_0(21)*occ_func_1_0(2)) + (occ_func_1_0(60)*occ_func_1_0(22)) + (occ_func_1_0(77)*occ_func_1_0(6)) + (occ_func_1_0(23)*occ_func_1_0(1)) + (occ_func_1_0(78)*occ_func_1_0(24)) + (occ_func_1_0(59)*occ_func_1_0(1)) + (occ_func_1_0(21)*occ_func_1_0(6)) + (occ_func_1_0(61)*occ_func_1_0(19)) + (occ_func_1_0(76)*occ_func_1_0(5)) + (occ_func_1_0(26)*occ_func_1_0(2)) + (occ_func_1_0(73)*occ_func_1_0(25)) + (occ_func_1_0(64)*occ_func_1_0(4)) + (occ_func_1_0(20)*occ_func_1_0(3)) + (occ_func_1_0(72)*occ_func_1_0(24)) + (occ_func_1_0(65)*occ_func_1_0(3)) + (occ_func_1_0(21)*occ_func_1_0(4)) + (occ_func_1_0(66)*occ_func_1_0(22)) + (occ_func_1_0(71)*occ_func_1_0(3)) + (occ_func_1_0(23)*occ_func_1_0(4)) + (occ_func_1_0(76)*occ_func_1_0(26)) + (occ_func_1_0(61)*occ_func_1_0(2)) + (occ_func_1_0(19)*occ_func_1_0(5)) + (occ_func_1_0(69)*occ_func_1_0(23)) + (occ_func_1_0(68)*occ_func_1_0(5)) + (occ_func_1_0(22)*occ_func_1_0(2)) + (occ_func_1_0(74)*occ_func_1_0(26)) + (occ_func_1_0(63)*occ_func_1_0(3)) + (occ_func_1_0(19)*occ_func_1_0(4)) + (occ_func_1_0(80)*occ_func_1_0(26)) + (occ_func_1_0(57)*occ_func_1_0(1)) + (occ_func_1_0(19)*occ_func_1_0(6)) + (occ_func_1_0(59)*occ_func_1_0(21)) + (occ_func_1_0(78)*occ_func_1_0(6)) + (occ_func_1_0(24)*occ_func_1_0(1)) + (occ_func_1_0(64)*occ_func_1_0(20)) + (occ_func_1_0(73)*occ_func_1_0(3)) + (occ_func_1_0(25)*occ_func_1_0(4)))/24.0;
  }

  /**** Basis functions for orbit 3, 14****
#Points: 3
MaxLength: 9.6513440  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
              -0.5000000   -1.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_14_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(61)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(77)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(62)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(65)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(59)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(69)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(64)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(73)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(79)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(60)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(63)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(78)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(66)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(72)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(57)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(75)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(70)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(68)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(74)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(71)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(80)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(76)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(67)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(58)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_14_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(61)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(76)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(77)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(60)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(62)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(75)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(65)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(72)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(59)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(78)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(69)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(68)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(64)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(73)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(73)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(64)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(79)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(58)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(60)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(77)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(63)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(74)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(78)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(59)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(66)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(71)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(72)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(65)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(57)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(80)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(75)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(62)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(70)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(67)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(68)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(69)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(74)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(63)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(71)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(66)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(80)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(57)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(76)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(61)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(67)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(70)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(58)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(79)*occ_func_1_0(17)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_14_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(61)) + (occ_func_1_0(74)*occ_func_1_0(12)) + (occ_func_1_0(76)*occ_func_1_0(13)) + (occ_func_1_0(71)*occ_func_1_0(77)) + (occ_func_1_0(66)*occ_func_1_0(17)) + (occ_func_1_0(60)*occ_func_1_0(8)) + (occ_func_1_0(58)*occ_func_1_0(62)) + (occ_func_1_0(79)*occ_func_1_0(15)) + (occ_func_1_0(75)*occ_func_1_0(10)) + (occ_func_1_0(67)*occ_func_1_0(65)) + (occ_func_1_0(70)*occ_func_1_0(11)) + (occ_func_1_0(72)*occ_func_1_0(14)) + (occ_func_1_0(65)*occ_func_1_0(59)) + (occ_func_1_0(72)*occ_func_1_0(9)) + (occ_func_1_0(78)*occ_func_1_0(16)) + (occ_func_1_0(77)*occ_func_1_0(69)) + (occ_func_1_0(60)*occ_func_1_0(7)) + (occ_func_1_0(68)*occ_func_1_0(18)) + (occ_func_1_0(62)*occ_func_1_0(64)) + (occ_func_1_0(75)*occ_func_1_0(14)) + (occ_func_1_0(73)*occ_func_1_0(11)) + (occ_func_1_0(79)*occ_func_1_0(73)) + (occ_func_1_0(58)*occ_func_1_0(8)) + (occ_func_1_0(64)*occ_func_1_0(17)) + (occ_func_1_0(75)*occ_func_1_0(79)) + (occ_func_1_0(62)*occ_func_1_0(15)) + (occ_func_1_0(58)*occ_func_1_0(10)) + (occ_func_1_0(68)*occ_func_1_0(60)) + (occ_func_1_0(69)*occ_func_1_0(7)) + (occ_func_1_0(77)*occ_func_1_0(18)) + (occ_func_1_0(57)*occ_func_1_0(63)) + (occ_func_1_0(80)*occ_func_1_0(16)) + (occ_func_1_0(74)*occ_func_1_0(9)) + (occ_func_1_0(70)*occ_func_1_0(78)) + (occ_func_1_0(67)*occ_func_1_0(18)) + (occ_func_1_0(59)*occ_func_1_0(7)) + (occ_func_1_0(60)*occ_func_1_0(66)) + (occ_func_1_0(77)*occ_func_1_0(17)) + (occ_func_1_0(71)*occ_func_1_0(8)) + (occ_func_1_0(78)*occ_func_1_0(72)) + (occ_func_1_0(59)*occ_func_1_0(9)) + (occ_func_1_0(65)*occ_func_1_0(16)) + (occ_func_1_0(61)*occ_func_1_0(57)) + (occ_func_1_0(76)*occ_func_1_0(10)) + (occ_func_1_0(80)*occ_func_1_0(15)) + (occ_func_1_0(73)*occ_func_1_0(75)) + (occ_func_1_0(64)*occ_func_1_0(14)) + (occ_func_1_0(62)*occ_func_1_0(11)) + (occ_func_1_0(72)*occ_func_1_0(70)) + (occ_func_1_0(65)*occ_func_1_0(11)) + (occ_func_1_0(67)*occ_func_1_0(14)) + (occ_func_1_0(66)*occ_func_1_0(68)) + (occ_func_1_0(71)*occ_func_1_0(13)) + (occ_func_1_0(69)*occ_func_1_0(12)) + (occ_func_1_0(76)*occ_func_1_0(74)) + (occ_func_1_0(61)*occ_func_1_0(12)) + (occ_func_1_0(63)*occ_func_1_0(13)) + (occ_func_1_0(69)*occ_func_1_0(71)) + (occ_func_1_0(68)*occ_func_1_0(13)) + (occ_func_1_0(66)*occ_func_1_0(12)) + (occ_func_1_0(74)*occ_func_1_0(80)) + (occ_func_1_0(63)*occ_func_1_0(16)) + (occ_func_1_0(57)*occ_func_1_0(9)) + (occ_func_1_0(80)*occ_func_1_0(76)) + (occ_func_1_0(57)*occ_func_1_0(10)) + (occ_func_1_0(61)*occ_func_1_0(15)) + (occ_func_1_0(59)*occ_func_1_0(67)) + (occ_func_1_0(78)*occ_func_1_0(18)) + (occ_func_1_0(70)*occ_func_1_0(7)) + (occ_func_1_0(64)*occ_func_1_0(58)) + (occ_func_1_0(73)*occ_func_1_0(8)) + (occ_func_1_0(79)*occ_func_1_0(17)))/24.0;
  }

  /**** Basis functions for orbit 3, 15****
#Points: 3
MaxLength: 9.6513440  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000    0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_15_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(29)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(30)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(30)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(31)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(28)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(30)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(32)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(27)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(30)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_15_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(29)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(30)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(29)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(30)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(27)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(32)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(31)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(28)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(29)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(30)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(32)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(27)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(28)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(31)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(32)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(27)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(31)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(28)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(31)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(28)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(27)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(32)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(28)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(31)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(27)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(32)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(32)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(27)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(28)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(31)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(29)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(30)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(30)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(29)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(30)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(29)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(31)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(28)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(28)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(31)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(30)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(7)) + (occ_func_1_0(29)*occ_func_1_0(18)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(32)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(27)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(27)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(32)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(30)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(29)*occ_func_1_0(7)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_15_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(29)) + (occ_func_1_0(74)*occ_func_1_0(18)) + (occ_func_1_0(30)*occ_func_1_0(7)) + (occ_func_1_0(71)*occ_func_1_0(29)) + (occ_func_1_0(66)*occ_func_1_0(10)) + (occ_func_1_0(30)*occ_func_1_0(15)) + (occ_func_1_0(58)*occ_func_1_0(27)) + (occ_func_1_0(79)*occ_func_1_0(13)) + (occ_func_1_0(32)*occ_func_1_0(12)) + (occ_func_1_0(67)*occ_func_1_0(31)) + (occ_func_1_0(70)*occ_func_1_0(17)) + (occ_func_1_0(28)*occ_func_1_0(8)) + (occ_func_1_0(65)*occ_func_1_0(29)) + (occ_func_1_0(72)*occ_func_1_0(15)) + (occ_func_1_0(30)*occ_func_1_0(10)) + (occ_func_1_0(77)*occ_func_1_0(32)) + (occ_func_1_0(60)*occ_func_1_0(14)) + (occ_func_1_0(27)*occ_func_1_0(11)) + (occ_func_1_0(62)*occ_func_1_0(28)) + (occ_func_1_0(75)*occ_func_1_0(16)) + (occ_func_1_0(31)*occ_func_1_0(9)) + (occ_func_1_0(79)*occ_func_1_0(32)) + (occ_func_1_0(58)*occ_func_1_0(12)) + (occ_func_1_0(27)*occ_func_1_0(13)) + (occ_func_1_0(75)*occ_func_1_0(31)) + (occ_func_1_0(62)*occ_func_1_0(9)) + (occ_func_1_0(28)*occ_func_1_0(16)) + (occ_func_1_0(68)*occ_func_1_0(31)) + (occ_func_1_0(69)*occ_func_1_0(16)) + (occ_func_1_0(28)*occ_func_1_0(9)) + (occ_func_1_0(57)*occ_func_1_0(27)) + (occ_func_1_0(80)*occ_func_1_0(14)) + (occ_func_1_0(32)*occ_func_1_0(11)) + (occ_func_1_0(70)*occ_func_1_0(28)) + (occ_func_1_0(67)*occ_func_1_0(8)) + (occ_func_1_0(31)*occ_func_1_0(17)) + (occ_func_1_0(60)*occ_func_1_0(27)) + (occ_func_1_0(77)*occ_func_1_0(11)) + (occ_func_1_0(32)*occ_func_1_0(14)) + (occ_func_1_0(78)*occ_func_1_0(32)) + (occ_func_1_0(59)*occ_func_1_0(13)) + (occ_func_1_0(27)*occ_func_1_0(12)) + (occ_func_1_0(61)*occ_func_1_0(28)) + (occ_func_1_0(76)*occ_func_1_0(17)) + (occ_func_1_0(31)*occ_func_1_0(8)) + (occ_func_1_0(73)*occ_func_1_0(29)) + (occ_func_1_0(64)*occ_func_1_0(7)) + (occ_func_1_0(30)*occ_func_1_0(18)) + (occ_func_1_0(72)*occ_func_1_0(30)) + (occ_func_1_0(65)*occ_func_1_0(10)) + (occ_func_1_0(29)*occ_func_1_0(15)) + (occ_func_1_0(66)*occ_func_1_0(30)) + (occ_func_1_0(71)*occ_func_1_0(15)) + (occ_func_1_0(29)*occ_func_1_0(10)) + (occ_func_1_0(76)*occ_func_1_0(31)) + (occ_func_1_0(61)*occ_func_1_0(8)) + (occ_func_1_0(28)*occ_func_1_0(17)) + (occ_func_1_0(69)*occ_func_1_0(28)) + (occ_func_1_0(68)*occ_func_1_0(9)) + (occ_func_1_0(31)*occ_func_1_0(16)) + (occ_func_1_0(74)*occ_func_1_0(30)) + (occ_func_1_0(63)*occ_func_1_0(7)) + (occ_func_1_0(29)*occ_func_1_0(18)) + (occ_func_1_0(80)*occ_func_1_0(32)) + (occ_func_1_0(57)*occ_func_1_0(11)) + (occ_func_1_0(27)*occ_func_1_0(14)) + (occ_func_1_0(59)*occ_func_1_0(27)) + (occ_func_1_0(78)*occ_func_1_0(12)) + (occ_func_1_0(32)*occ_func_1_0(13)) + (occ_func_1_0(64)*occ_func_1_0(30)) + (occ_func_1_0(73)*occ_func_1_0(18)) + (occ_func_1_0(29)*occ_func_1_0(7)))/24.0;
  }

  /**** Basis functions for orbit 3, 16****
#Points: 3
MaxLength: 9.6513440  MinLength: 5.5722061
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000   -0.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_16_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(16)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(15)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(17)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(12)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(13)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(11)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(14)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(8)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(15)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(9)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(7)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(10)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(16)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_16_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(11)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(14)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(16)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(9)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(7)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(18)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(13)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(12)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(8)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(17)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(15)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(10)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(12)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(13)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(10)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(11)) + (occ_func_1_0(15)*occ_func_1_0(14)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(8)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(17)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(17)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(8)*occ_func_1_0(15)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(7)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(18)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(13)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(12)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(12)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(16)) + (occ_func_1_0(13)*occ_func_1_0(9)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(13)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(12)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(11)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(9)) + (occ_func_1_0(14)*occ_func_1_0(16)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(11)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(14)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(14)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(15)) + (occ_func_1_0(11)*occ_func_1_0(10)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(8)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(14)) + (occ_func_1_0(17)*occ_func_1_0(11)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(7)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(17)) + (occ_func_1_0(18)*occ_func_1_0(8)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(15)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(8)) + (occ_func_1_0(10)*occ_func_1_0(17)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(9)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(18)) + (occ_func_1_0(16)*occ_func_1_0(7)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(7)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(13)) + (occ_func_1_0(18)*occ_func_1_0(12)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(10)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(12)) + (occ_func_1_0(15)*occ_func_1_0(13)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(16)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(10)) + (occ_func_1_0(9)*occ_func_1_0(15)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_16_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(11)) + (occ_func_1_0(74)*occ_func_1_0(17)) + (occ_func_1_0(14)*occ_func_1_0(8)) + (occ_func_1_0(71)*occ_func_1_0(16)) + (occ_func_1_0(66)*occ_func_1_0(14)) + (occ_func_1_0(9)*occ_func_1_0(11)) + (occ_func_1_0(58)*occ_func_1_0(7)) + (occ_func_1_0(79)*occ_func_1_0(16)) + (occ_func_1_0(18)*occ_func_1_0(9)) + (occ_func_1_0(67)*occ_func_1_0(13)) + (occ_func_1_0(70)*occ_func_1_0(15)) + (occ_func_1_0(12)*occ_func_1_0(10)) + (occ_func_1_0(65)*occ_func_1_0(8)) + (occ_func_1_0(72)*occ_func_1_0(12)) + (occ_func_1_0(17)*occ_func_1_0(13)) + (occ_func_1_0(77)*occ_func_1_0(15)) + (occ_func_1_0(60)*occ_func_1_0(9)) + (occ_func_1_0(10)*occ_func_1_0(16)) + (occ_func_1_0(62)*occ_func_1_0(12)) + (occ_func_1_0(75)*occ_func_1_0(18)) + (occ_func_1_0(13)*occ_func_1_0(7)) + (occ_func_1_0(68)*occ_func_1_0(10)) + (occ_func_1_0(69)*occ_func_1_0(11)) + (occ_func_1_0(15)*occ_func_1_0(14)) + (occ_func_1_0(57)*occ_func_1_0(8)) + (occ_func_1_0(80)*occ_func_1_0(18)) + (occ_func_1_0(17)*occ_func_1_0(7)) + (occ_func_1_0(78)*occ_func_1_0(17)) + (occ_func_1_0(59)*occ_func_1_0(10)) + (occ_func_1_0(8)*occ_func_1_0(15)) + (occ_func_1_0(61)*occ_func_1_0(7)) + (occ_func_1_0(76)*occ_func_1_0(14)) + (occ_func_1_0(18)*occ_func_1_0(11)) + (occ_func_1_0(73)*occ_func_1_0(13)) + (occ_func_1_0(64)*occ_func_1_0(9)) + (occ_func_1_0(12)*occ_func_1_0(16)) + (occ_func_1_0(64)*occ_func_1_0(12)) + (occ_func_1_0(73)*occ_func_1_0(16)) + (occ_func_1_0(13)*occ_func_1_0(9)) + (occ_func_1_0(65)*occ_func_1_0(13)) + (occ_func_1_0(72)*occ_func_1_0(17)) + (occ_func_1_0(12)*occ_func_1_0(8)) + (occ_func_1_0(71)*occ_func_1_0(11)) + (occ_func_1_0(66)*occ_func_1_0(9)) + (occ_func_1_0(14)*occ_func_1_0(16)) + (occ_func_1_0(61)*occ_func_1_0(11)) + (occ_func_1_0(76)*occ_func_1_0(18)) + (occ_func_1_0(14)*occ_func_1_0(7)) + (occ_func_1_0(68)*occ_func_1_0(14)) + (occ_func_1_0(69)*occ_func_1_0(15)) + (occ_func_1_0(11)*occ_func_1_0(10)) + (occ_func_1_0(63)*occ_func_1_0(8)) + (occ_func_1_0(74)*occ_func_1_0(14)) + (occ_func_1_0(17)*occ_func_1_0(11)) + (occ_func_1_0(57)*occ_func_1_0(7)) + (occ_func_1_0(80)*occ_func_1_0(17)) + (occ_func_1_0(18)*occ_func_1_0(8)) + (occ_func_1_0(78)*occ_func_1_0(15)) + (occ_func_1_0(59)*occ_func_1_0(8)) + (occ_func_1_0(10)*occ_func_1_0(17)) + (occ_func_1_0(58)*occ_func_1_0(9)) + (occ_func_1_0(79)*occ_func_1_0(18)) + (occ_func_1_0(16)*occ_func_1_0(7)) + (occ_func_1_0(62)*occ_func_1_0(7)) + (occ_func_1_0(75)*occ_func_1_0(13)) + (occ_func_1_0(18)*occ_func_1_0(12)) + (occ_func_1_0(67)*occ_func_1_0(10)) + (occ_func_1_0(70)*occ_func_1_0(12)) + (occ_func_1_0(15)*occ_func_1_0(13)) + (occ_func_1_0(77)*occ_func_1_0(16)) + (occ_func_1_0(60)*occ_func_1_0(10)) + (occ_func_1_0(9)*occ_func_1_0(15)))/24.0;
  }

  /**** Basis functions for orbit 3, 17****
#Points: 3
MaxLength: 9.6513440  MinLength: 6.8245308
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
               0.5000000   -1.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_17_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(44)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(39)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(42)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(48)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(46)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(42)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(46)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(39)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(44)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(48)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_17_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(41)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(48)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(54)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(35)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(37)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(52)*occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(45)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(44)*occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(34)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(55)*occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(49)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(40)*occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(44)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(45)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(50)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(39)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(56)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(33)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(36)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(53)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(38)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(25)) + (occ_func_1_0(51)*occ_func_1_0(20)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(53)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(25)) + (occ_func_1_0(36)*occ_func_1_0(20)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(39)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(50)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(51)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(38)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(33)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(56)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(47)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(42)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(42)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(47)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(48)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(25)) + (occ_func_1_0(41)*occ_func_1_0(20)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(46)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(43)*occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(43)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(46)*occ_func_1_0(24)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(55)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(34)*occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(52)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(37)*occ_func_1_0(24)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(40)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(49)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(35)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(54)*occ_func_1_0(24)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(42)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(47)*occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(47)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(42)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(41)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(48)*occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(43)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(25)) + (occ_func_1_0(46)*occ_func_1_0(20)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(46)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(43)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(34)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(55)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(37)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(52)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(49)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(40)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(54)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(24)) + (occ_func_1_0(35)*occ_func_1_0(21)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(39)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(50)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(33)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(56)*occ_func_1_0(24)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(53)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(26)) + (occ_func_1_0(36)*occ_func_1_0(19)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(51)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(38)*occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(36)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(20)) + (occ_func_1_0(53)*occ_func_1_0(25)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(50)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(39)*occ_func_1_0(24)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(38)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(51)*occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(56)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(23)) + (occ_func_1_0(33)*occ_func_1_0(22)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(35)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(54)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(52)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(37)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(44)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(22)) + (occ_func_1_0(45)*occ_func_1_0(23)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(55)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(25)) + (occ_func_1_0(34)*occ_func_1_0(20)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(40)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(25)) + (occ_func_1_0(49)*occ_func_1_0(20)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(45)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(19)) + (occ_func_1_0(44)*occ_func_1_0(26)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(48)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(21)) + (occ_func_1_0(41)*occ_func_1_0(24)*occ_func_1_0(0)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_17_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(41)) + (occ_func_1_0(74)*occ_func_1_0(24)) + (occ_func_1_0(48)*occ_func_1_0(21)) + (occ_func_1_0(71)*occ_func_1_0(54)) + (occ_func_1_0(66)*occ_func_1_0(26)) + (occ_func_1_0(35)*occ_func_1_0(19)) + (occ_func_1_0(58)*occ_func_1_0(37)) + (occ_func_1_0(79)*occ_func_1_0(23)) + (occ_func_1_0(52)*occ_func_1_0(22)) + (occ_func_1_0(67)*occ_func_1_0(45)) + (occ_func_1_0(70)*occ_func_1_0(23)) + (occ_func_1_0(44)*occ_func_1_0(22)) + (occ_func_1_0(65)*occ_func_1_0(34)) + (occ_func_1_0(72)*occ_func_1_0(20)) + (occ_func_1_0(55)*occ_func_1_0(25)) + (occ_func_1_0(77)*occ_func_1_0(49)) + (occ_func_1_0(60)*occ_func_1_0(20)) + (occ_func_1_0(40)*occ_func_1_0(25)) + (occ_func_1_0(62)*occ_func_1_0(44)) + (occ_func_1_0(75)*occ_func_1_0(26)) + (occ_func_1_0(45)*occ_func_1_0(19)) + (occ_func_1_0(79)*occ_func_1_0(50)) + (occ_func_1_0(58)*occ_func_1_0(19)) + (occ_func_1_0(39)*occ_func_1_0(26)) + (occ_func_1_0(75)*occ_func_1_0(56)) + (occ_func_1_0(62)*occ_func_1_0(24)) + (occ_func_1_0(33)*occ_func_1_0(21)) + (occ_func_1_0(68)*occ_func_1_0(36)) + (occ_func_1_0(69)*occ_func_1_0(19)) + (occ_func_1_0(53)*occ_func_1_0(26)) + (occ_func_1_0(57)*occ_func_1_0(38)) + (occ_func_1_0(80)*occ_func_1_0(25)) + (occ_func_1_0(51)*occ_func_1_0(20)) + (occ_func_1_0(70)*occ_func_1_0(53)) + (occ_func_1_0(67)*occ_func_1_0(25)) + (occ_func_1_0(36)*occ_func_1_0(20)) + (occ_func_1_0(60)*occ_func_1_0(39)) + (occ_func_1_0(77)*occ_func_1_0(24)) + (occ_func_1_0(50)*occ_func_1_0(21)) + (occ_func_1_0(78)*occ_func_1_0(51)) + (occ_func_1_0(59)*occ_func_1_0(22)) + (occ_func_1_0(38)*occ_func_1_0(23)) + (occ_func_1_0(61)*occ_func_1_0(33)) + (occ_func_1_0(76)*occ_func_1_0(22)) + (occ_func_1_0(56)*occ_func_1_0(23)) + (occ_func_1_0(73)*occ_func_1_0(47)) + (occ_func_1_0(64)*occ_func_1_0(22)) + (occ_func_1_0(42)*occ_func_1_0(23)) + (occ_func_1_0(72)*occ_func_1_0(42)) + (occ_func_1_0(65)*occ_func_1_0(19)) + (occ_func_1_0(47)*occ_func_1_0(26)) + (occ_func_1_0(66)*occ_func_1_0(48)) + (occ_func_1_0(71)*occ_func_1_0(25)) + (occ_func_1_0(41)*occ_func_1_0(20)) + (occ_func_1_0(76)*occ_func_1_0(46)) + (occ_func_1_0(61)*occ_func_1_0(20)) + (occ_func_1_0(43)*occ_func_1_0(25)) + (occ_func_1_0(69)*occ_func_1_0(43)) + (occ_func_1_0(68)*occ_func_1_0(21)) + (occ_func_1_0(46)*occ_func_1_0(24)) + (occ_func_1_0(74)*occ_func_1_0(55)) + (occ_func_1_0(63)*occ_func_1_0(23)) + (occ_func_1_0(34)*occ_func_1_0(22)) + (occ_func_1_0(80)*occ_func_1_0(52)) + (occ_func_1_0(57)*occ_func_1_0(21)) + (occ_func_1_0(37)*occ_func_1_0(24)) + (occ_func_1_0(59)*occ_func_1_0(40)) + (occ_func_1_0(78)*occ_func_1_0(26)) + (occ_func_1_0(49)*occ_func_1_0(19)) + (occ_func_1_0(64)*occ_func_1_0(35)) + (occ_func_1_0(73)*occ_func_1_0(21)) + (occ_func_1_0(54)*occ_func_1_0(24)) + (occ_func_1_0(64)*occ_func_1_0(42)) + (occ_func_1_0(73)*occ_func_1_0(23)) + (occ_func_1_0(47)*occ_func_1_0(22)) + (occ_func_1_0(65)*occ_func_1_0(47)) + (occ_func_1_0(72)*occ_func_1_0(26)) + (occ_func_1_0(42)*occ_func_1_0(19)) + (occ_func_1_0(71)*occ_func_1_0(41)) + (occ_func_1_0(66)*occ_func_1_0(20)) + (occ_func_1_0(48)*occ_func_1_0(25)) + (occ_func_1_0(61)*occ_func_1_0(43)) + (occ_func_1_0(76)*occ_func_1_0(25)) + (occ_func_1_0(46)*occ_func_1_0(20)) + (occ_func_1_0(68)*occ_func_1_0(46)) + (occ_func_1_0(69)*occ_func_1_0(24)) + (occ_func_1_0(43)*occ_func_1_0(21)) + (occ_func_1_0(63)*occ_func_1_0(34)) + (occ_func_1_0(74)*occ_func_1_0(22)) + (occ_func_1_0(55)*occ_func_1_0(23)) + (occ_func_1_0(57)*occ_func_1_0(37)) + (occ_func_1_0(80)*occ_func_1_0(24)) + (occ_func_1_0(52)*occ_func_1_0(21)) + (occ_func_1_0(78)*occ_func_1_0(49)) + (occ_func_1_0(59)*occ_func_1_0(19)) + (occ_func_1_0(40)*occ_func_1_0(26)) + (occ_func_1_0(73)*occ_func_1_0(54)) + (occ_func_1_0(64)*occ_func_1_0(24)) + (occ_func_1_0(35)*occ_func_1_0(21)) + (occ_func_1_0(58)*occ_func_1_0(39)) + (occ_func_1_0(79)*occ_func_1_0(26)) + (occ_func_1_0(50)*occ_func_1_0(19)) + (occ_func_1_0(62)*occ_func_1_0(33)) + (occ_func_1_0(75)*occ_func_1_0(21)) + (occ_func_1_0(56)*occ_func_1_0(24)) + (occ_func_1_0(69)*occ_func_1_0(53)) + (occ_func_1_0(68)*occ_func_1_0(26)) + (occ_func_1_0(36)*occ_func_1_0(19)) + (occ_func_1_0(80)*occ_func_1_0(51)) + (occ_func_1_0(57)*occ_func_1_0(20)) + (occ_func_1_0(38)*occ_func_1_0(25)) + (occ_func_1_0(67)*occ_func_1_0(36)) + (occ_func_1_0(70)*occ_func_1_0(20)) + (occ_func_1_0(53)*occ_func_1_0(25)) + (occ_func_1_0(77)*occ_func_1_0(50)) + (occ_func_1_0(60)*occ_func_1_0(21)) + (occ_func_1_0(39)*occ_func_1_0(24)) + (occ_func_1_0(59)*occ_func_1_0(38)) + (occ_func_1_0(78)*occ_func_1_0(23)) + (occ_func_1_0(51)*occ_func_1_0(22)) + (occ_func_1_0(76)*occ_func_1_0(56)) + (occ_func_1_0(61)*occ_func_1_0(23)) + (occ_func_1_0(33)*occ_func_1_0(22)) + (occ_func_1_0(66)*occ_func_1_0(35)) + (occ_func_1_0(71)*occ_func_1_0(19)) + (occ_func_1_0(54)*occ_func_1_0(26)) + (occ_func_1_0(79)*occ_func_1_0(52)) + (occ_func_1_0(58)*occ_func_1_0(22)) + (occ_func_1_0(37)*occ_func_1_0(23)) + (occ_func_1_0(70)*occ_func_1_0(44)) + (occ_func_1_0(67)*occ_func_1_0(22)) + (occ_func_1_0(45)*occ_func_1_0(23)) + (occ_func_1_0(72)*occ_func_1_0(55)) + (occ_func_1_0(65)*occ_func_1_0(25)) + (occ_func_1_0(34)*occ_func_1_0(20)) + (occ_func_1_0(60)*occ_func_1_0(40)) + (occ_func_1_0(77)*occ_func_1_0(25)) + (occ_func_1_0(49)*occ_func_1_0(20)) + (occ_func_1_0(75)*occ_func_1_0(45)) + (occ_func_1_0(62)*occ_func_1_0(19)) + (occ_func_1_0(44)*occ_func_1_0(26)) + (occ_func_1_0(74)*occ_func_1_0(48)) + (occ_func_1_0(63)*occ_func_1_0(21)) + (occ_func_1_0(41)*occ_func_1_0(24)))/48.0;
  }

  /**** Basis functions for orbit 3, 18****
#Points: 3
MaxLength: 9.6513440  MinLength: 7.8802894
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
               1.5000000   -0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_18_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(71)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(73)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(57)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(75)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(63)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(78)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(70)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(77)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(76)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(67)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(59)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(69)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(58)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(80)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(62)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(65)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(64)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(74)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(68)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(61)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(72)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(79)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(60)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(66)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_18_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(71)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(32)) + (occ_func_1_0(66)*occ_func_1_0(27)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(73)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(64)*occ_func_1_0(28)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(57)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(29)) + (occ_func_1_0(80)*occ_func_1_0(30)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(75)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(32)) + (occ_func_1_0(62)*occ_func_1_0(27)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(63)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(28)) + (occ_func_1_0(74)*occ_func_1_0(31)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(78)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(30)) + (occ_func_1_0(59)*occ_func_1_0(29)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(70)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(32)) + (occ_func_1_0(67)*occ_func_1_0(27)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(77)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(28)) + (occ_func_1_0(60)*occ_func_1_0(31)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(76)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(30)) + (occ_func_1_0(61)*occ_func_1_0(29)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(67)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(29)) + (occ_func_1_0(70)*occ_func_1_0(30)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(59)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(78)*occ_func_1_0(28)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(69)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(29)) + (occ_func_1_0(68)*occ_func_1_0(30)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(58)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(28)) + (occ_func_1_0(79)*occ_func_1_0(31)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(80)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(57)*occ_func_1_0(28)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(62)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(30)) + (occ_func_1_0(75)*occ_func_1_0(29)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(65)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(72)*occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(64)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(73)*occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(74)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(32)) + (occ_func_1_0(63)*occ_func_1_0(27)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(68)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(69)*occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(61)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(27)) + (occ_func_1_0(76)*occ_func_1_0(32)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(72)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(28)) + (occ_func_1_0(65)*occ_func_1_0(31)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(79)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(29)) + (occ_func_1_0(58)*occ_func_1_0(30)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(60)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(30)) + (occ_func_1_0(77)*occ_func_1_0(29)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(66)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(31)) + (occ_func_1_0(71)*occ_func_1_0(28)*occ_func_1_0(0)))/24.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_18_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(71)) + (occ_func_1_0(74)*occ_func_1_0(32)) + (occ_func_1_0(66)*occ_func_1_0(27)) + (occ_func_1_0(71)*occ_func_1_0(73)) + (occ_func_1_0(66)*occ_func_1_0(31)) + (occ_func_1_0(64)*occ_func_1_0(28)) + (occ_func_1_0(58)*occ_func_1_0(57)) + (occ_func_1_0(79)*occ_func_1_0(29)) + (occ_func_1_0(80)*occ_func_1_0(30)) + (occ_func_1_0(67)*occ_func_1_0(75)) + (occ_func_1_0(70)*occ_func_1_0(32)) + (occ_func_1_0(62)*occ_func_1_0(27)) + (occ_func_1_0(65)*occ_func_1_0(63)) + (occ_func_1_0(72)*occ_func_1_0(28)) + (occ_func_1_0(74)*occ_func_1_0(31)) + (occ_func_1_0(77)*occ_func_1_0(78)) + (occ_func_1_0(60)*occ_func_1_0(30)) + (occ_func_1_0(59)*occ_func_1_0(29)) + (occ_func_1_0(62)*occ_func_1_0(70)) + (occ_func_1_0(75)*occ_func_1_0(32)) + (occ_func_1_0(67)*occ_func_1_0(27)) + (occ_func_1_0(79)*occ_func_1_0(77)) + (occ_func_1_0(58)*occ_func_1_0(28)) + (occ_func_1_0(60)*occ_func_1_0(31)) + (occ_func_1_0(75)*occ_func_1_0(76)) + (occ_func_1_0(62)*occ_func_1_0(30)) + (occ_func_1_0(61)*occ_func_1_0(29)) + (occ_func_1_0(68)*occ_func_1_0(67)) + (occ_func_1_0(69)*occ_func_1_0(29)) + (occ_func_1_0(70)*occ_func_1_0(30)) + (occ_func_1_0(57)*occ_func_1_0(59)) + (occ_func_1_0(80)*occ_func_1_0(31)) + (occ_func_1_0(78)*occ_func_1_0(28)) + (occ_func_1_0(70)*occ_func_1_0(69)) + (occ_func_1_0(67)*occ_func_1_0(29)) + (occ_func_1_0(68)*occ_func_1_0(30)) + (occ_func_1_0(60)*occ_func_1_0(58)) + (occ_func_1_0(77)*occ_func_1_0(28)) + (occ_func_1_0(79)*occ_func_1_0(31)) + (occ_func_1_0(78)*occ_func_1_0(80)) + (occ_func_1_0(59)*occ_func_1_0(31)) + (occ_func_1_0(57)*occ_func_1_0(28)) + (occ_func_1_0(61)*occ_func_1_0(62)) + (occ_func_1_0(76)*occ_func_1_0(30)) + (occ_func_1_0(75)*occ_func_1_0(29)) + (occ_func_1_0(73)*occ_func_1_0(65)) + (occ_func_1_0(64)*occ_func_1_0(27)) + (occ_func_1_0(72)*occ_func_1_0(32)) + (occ_func_1_0(72)*occ_func_1_0(64)) + (occ_func_1_0(65)*occ_func_1_0(27)) + (occ_func_1_0(73)*occ_func_1_0(32)) + (occ_func_1_0(66)*occ_func_1_0(74)) + (occ_func_1_0(71)*occ_func_1_0(32)) + (occ_func_1_0(63)*occ_func_1_0(27)) + (occ_func_1_0(76)*occ_func_1_0(68)) + (occ_func_1_0(61)*occ_func_1_0(27)) + (occ_func_1_0(69)*occ_func_1_0(32)) + (occ_func_1_0(69)*occ_func_1_0(61)) + (occ_func_1_0(68)*occ_func_1_0(27)) + (occ_func_1_0(76)*occ_func_1_0(32)) + (occ_func_1_0(74)*occ_func_1_0(72)) + (occ_func_1_0(63)*occ_func_1_0(28)) + (occ_func_1_0(65)*occ_func_1_0(31)) + (occ_func_1_0(80)*occ_func_1_0(79)) + (occ_func_1_0(57)*occ_func_1_0(29)) + (occ_func_1_0(58)*occ_func_1_0(30)) + (occ_func_1_0(59)*occ_func_1_0(60)) + (occ_func_1_0(78)*occ_func_1_0(30)) + (occ_func_1_0(77)*occ_func_1_0(29)) + (occ_func_1_0(64)*occ_func_1_0(66)) + (occ_func_1_0(73)*occ_func_1_0(31)) + (occ_func_1_0(71)*occ_func_1_0(28)))/24.0;
  }

  /**** Basis functions for orbit 3, 19****
#Points: 3
MaxLength: 9.6513440  MinLength: 8.8104314
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
               1.5000000    0.5000000   -1.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_19_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(48)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(42)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(39)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(44)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(46)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(51)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(50)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(38)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(49)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(52)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(45)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(35)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(54)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(43)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(36)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(41)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(42)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(53)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(48)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(56)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(33)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(47)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(44)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(55)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(37)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(46)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(34)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(40)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(39)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_19_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(50)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(56)) + (occ_func_1_0(39)*occ_func_1_0(33)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(45)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(40)) + (occ_func_1_0(44)*occ_func_1_0(49)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(34)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(45)) + (occ_func_1_0(55)*occ_func_1_0(44)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(52)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(55)) + (occ_func_1_0(37)*occ_func_1_0(34)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(43)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(49)) + (occ_func_1_0(46)*occ_func_1_0(40)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(55)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(46)) + (occ_func_1_0(34)*occ_func_1_0(43)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(49)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(54)) + (occ_func_1_0(40)*occ_func_1_0(35)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(53)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(42)) + (occ_func_1_0(36)*occ_func_1_0(47)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(48)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(39)) + (occ_func_1_0(41)*occ_func_1_0(50)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(47)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(50)) + (occ_func_1_0(42)*occ_func_1_0(39)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(36)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(48)) + (occ_func_1_0(53)*occ_func_1_0(41)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(41)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(38)) + (occ_func_1_0(48)*occ_func_1_0(51)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(33)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(41)) + (occ_func_1_0(56)*occ_func_1_0(48)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(56)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(47)) + (occ_func_1_0(33)*occ_func_1_0(42)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(42)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(51)) + (occ_func_1_0(47)*occ_func_1_0(38)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(38)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(33)) + (occ_func_1_0(51)*occ_func_1_0(56)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(39)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(36)) + (occ_func_1_0(50)*occ_func_1_0(53)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(51)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(53)) + (occ_func_1_0(38)*occ_func_1_0(36)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(40)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(34)) + (occ_func_1_0(49)*occ_func_1_0(55)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(37)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(35)) + (occ_func_1_0(52)*occ_func_1_0(54)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(44)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(37)) + (occ_func_1_0(45)*occ_func_1_0(52)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(54)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(43)) + (occ_func_1_0(35)*occ_func_1_0(46)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(35)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(44)) + (occ_func_1_0(54)*occ_func_1_0(45)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(46)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(52)) + (occ_func_1_0(43)*occ_func_1_0(37)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(51)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(56)) + (occ_func_1_0(38)*occ_func_1_0(33)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(50)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(53)) + (occ_func_1_0(39)*occ_func_1_0(36)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(38)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(36)) + (occ_func_1_0(51)*occ_func_1_0(53)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(61)*occ_func_1_0(49)) + (occ_func_1_0(76)*occ_func_1_0(0)*occ_func_1_0(55)) + (occ_func_1_0(40)*occ_func_1_0(34)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(52)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(54)) + (occ_func_1_0(37)*occ_func_1_0(35)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(45)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(52)) + (occ_func_1_0(44)*occ_func_1_0(37)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(57)*occ_func_1_0(35)) + (occ_func_1_0(80)*occ_func_1_0(0)*occ_func_1_0(46)) + (occ_func_1_0(54)*occ_func_1_0(43)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(78)*occ_func_1_0(54)) + (occ_func_1_0(59)*occ_func_1_0(0)*occ_func_1_0(45)) + (occ_func_1_0(35)*occ_func_1_0(44)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(73)*occ_func_1_0(43)) + (occ_func_1_0(64)*occ_func_1_0(0)*occ_func_1_0(37)) + (occ_func_1_0(46)*occ_func_1_0(52)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(58)*occ_func_1_0(36)) + (occ_func_1_0(79)*occ_func_1_0(0)*occ_func_1_0(47)) + (occ_func_1_0(53)*occ_func_1_0(42)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(62)*occ_func_1_0(41)) + (occ_func_1_0(75)*occ_func_1_0(0)*occ_func_1_0(50)) + (occ_func_1_0(48)*occ_func_1_0(39)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(69)*occ_func_1_0(42)) + (occ_func_1_0(68)*occ_func_1_0(0)*occ_func_1_0(39)) + (occ_func_1_0(47)*occ_func_1_0(50)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(80)*occ_func_1_0(53)) + (occ_func_1_0(57)*occ_func_1_0(0)*occ_func_1_0(41)) + (occ_func_1_0(36)*occ_func_1_0(48)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(48)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(51)) + (occ_func_1_0(41)*occ_func_1_0(38)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(77)*occ_func_1_0(56)) + (occ_func_1_0(60)*occ_func_1_0(0)*occ_func_1_0(48)) + (occ_func_1_0(33)*occ_func_1_0(41)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(59)*occ_func_1_0(33)) + (occ_func_1_0(78)*occ_func_1_0(0)*occ_func_1_0(42)) + (occ_func_1_0(56)*occ_func_1_0(47)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(76)*occ_func_1_0(47)) + (occ_func_1_0(61)*occ_func_1_0(0)*occ_func_1_0(38)) + (occ_func_1_0(42)*occ_func_1_0(51)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(66)*occ_func_1_0(44)) + (occ_func_1_0(71)*occ_func_1_0(0)*occ_func_1_0(49)) + (occ_func_1_0(45)*occ_func_1_0(40)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(79)*occ_func_1_0(55)) + (occ_func_1_0(58)*occ_func_1_0(0)*occ_func_1_0(44)) + (occ_func_1_0(34)*occ_func_1_0(45)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(70)*occ_func_1_0(37)) + (occ_func_1_0(67)*occ_func_1_0(0)*occ_func_1_0(34)) + (occ_func_1_0(52)*occ_func_1_0(55)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(72)*occ_func_1_0(46)) + (occ_func_1_0(65)*occ_func_1_0(0)*occ_func_1_0(40)) + (occ_func_1_0(43)*occ_func_1_0(49)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(60)*occ_func_1_0(34)) + (occ_func_1_0(77)*occ_func_1_0(0)*occ_func_1_0(43)) + (occ_func_1_0(55)*occ_func_1_0(46)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(75)*occ_func_1_0(40)) + (occ_func_1_0(62)*occ_func_1_0(0)*occ_func_1_0(35)) + (occ_func_1_0(49)*occ_func_1_0(54)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(74)*occ_func_1_0(39)) + (occ_func_1_0(63)*occ_func_1_0(0)*occ_func_1_0(33)) + (occ_func_1_0(50)*occ_func_1_0(56)*occ_func_1_0(0)))/48.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_19_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(50)) + (occ_func_1_0(74)*occ_func_1_0(56)) + (occ_func_1_0(39)*occ_func_1_0(33)) + (occ_func_1_0(71)*occ_func_1_0(45)) + (occ_func_1_0(66)*occ_func_1_0(40)) + (occ_func_1_0(44)*occ_func_1_0(49)) + (occ_func_1_0(58)*occ_func_1_0(34)) + (occ_func_1_0(79)*occ_func_1_0(45)) + (occ_func_1_0(55)*occ_func_1_0(44)) + (occ_func_1_0(67)*occ_func_1_0(52)) + (occ_func_1_0(70)*occ_func_1_0(55)) + (occ_func_1_0(37)*occ_func_1_0(34)) + (occ_func_1_0(65)*occ_func_1_0(43)) + (occ_func_1_0(72)*occ_func_1_0(49)) + (occ_func_1_0(46)*occ_func_1_0(40)) + (occ_func_1_0(77)*occ_func_1_0(55)) + (occ_func_1_0(60)*occ_func_1_0(46)) + (occ_func_1_0(34)*occ_func_1_0(43)) + (occ_func_1_0(62)*occ_func_1_0(49)) + (occ_func_1_0(75)*occ_func_1_0(54)) + (occ_func_1_0(40)*occ_func_1_0(35)) + (occ_func_1_0(79)*occ_func_1_0(53)) + (occ_func_1_0(58)*occ_func_1_0(42)) + (occ_func_1_0(36)*occ_func_1_0(47)) + (occ_func_1_0(75)*occ_func_1_0(48)) + (occ_func_1_0(62)*occ_func_1_0(39)) + (occ_func_1_0(41)*occ_func_1_0(50)) + (occ_func_1_0(68)*occ_func_1_0(47)) + (occ_func_1_0(69)*occ_func_1_0(50)) + (occ_func_1_0(42)*occ_func_1_0(39)) + (occ_func_1_0(57)*occ_func_1_0(36)) + (occ_func_1_0(80)*occ_func_1_0(48)) + (occ_func_1_0(53)*occ_func_1_0(41)) + (occ_func_1_0(70)*occ_func_1_0(41)) + (occ_func_1_0(67)*occ_func_1_0(38)) + (occ_func_1_0(48)*occ_func_1_0(51)) + (occ_func_1_0(60)*occ_func_1_0(33)) + (occ_func_1_0(77)*occ_func_1_0(41)) + (occ_func_1_0(56)*occ_func_1_0(48)) + (occ_func_1_0(78)*occ_func_1_0(56)) + (occ_func_1_0(59)*occ_func_1_0(47)) + (occ_func_1_0(33)*occ_func_1_0(42)) + (occ_func_1_0(61)*occ_func_1_0(42)) + (occ_func_1_0(76)*occ_func_1_0(51)) + (occ_func_1_0(47)*occ_func_1_0(38)) + (occ_func_1_0(73)*occ_func_1_0(38)) + (occ_func_1_0(64)*occ_func_1_0(33)) + (occ_func_1_0(51)*occ_func_1_0(56)) + (occ_func_1_0(72)*occ_func_1_0(39)) + (occ_func_1_0(65)*occ_func_1_0(36)) + (occ_func_1_0(50)*occ_func_1_0(53)) + (occ_func_1_0(66)*occ_func_1_0(51)) + (occ_func_1_0(71)*occ_func_1_0(53)) + (occ_func_1_0(38)*occ_func_1_0(36)) + (occ_func_1_0(76)*occ_func_1_0(40)) + (occ_func_1_0(61)*occ_func_1_0(34)) + (occ_func_1_0(49)*occ_func_1_0(55)) + (occ_func_1_0(69)*occ_func_1_0(37)) + (occ_func_1_0(68)*occ_func_1_0(35)) + (occ_func_1_0(52)*occ_func_1_0(54)) + (occ_func_1_0(74)*occ_func_1_0(44)) + (occ_func_1_0(63)*occ_func_1_0(37)) + (occ_func_1_0(45)*occ_func_1_0(52)) + (occ_func_1_0(80)*occ_func_1_0(54)) + (occ_func_1_0(57)*occ_func_1_0(43)) + (occ_func_1_0(35)*occ_func_1_0(46)) + (occ_func_1_0(59)*occ_func_1_0(35)) + (occ_func_1_0(78)*occ_func_1_0(44)) + (occ_func_1_0(54)*occ_func_1_0(45)) + (occ_func_1_0(64)*occ_func_1_0(46)) + (occ_func_1_0(73)*occ_func_1_0(52)) + (occ_func_1_0(43)*occ_func_1_0(37)) + (occ_func_1_0(64)*occ_func_1_0(51)) + (occ_func_1_0(73)*occ_func_1_0(56)) + (occ_func_1_0(38)*occ_func_1_0(33)) + (occ_func_1_0(65)*occ_func_1_0(50)) + (occ_func_1_0(72)*occ_func_1_0(53)) + (occ_func_1_0(39)*occ_func_1_0(36)) + (occ_func_1_0(71)*occ_func_1_0(38)) + (occ_func_1_0(66)*occ_func_1_0(36)) + (occ_func_1_0(51)*occ_func_1_0(53)) + (occ_func_1_0(61)*occ_func_1_0(49)) + (occ_func_1_0(76)*occ_func_1_0(55)) + (occ_func_1_0(40)*occ_func_1_0(34)) + (occ_func_1_0(68)*occ_func_1_0(52)) + (occ_func_1_0(69)*occ_func_1_0(54)) + (occ_func_1_0(37)*occ_func_1_0(35)) + (occ_func_1_0(63)*occ_func_1_0(45)) + (occ_func_1_0(74)*occ_func_1_0(52)) + (occ_func_1_0(44)*occ_func_1_0(37)) + (occ_func_1_0(57)*occ_func_1_0(35)) + (occ_func_1_0(80)*occ_func_1_0(46)) + (occ_func_1_0(54)*occ_func_1_0(43)) + (occ_func_1_0(78)*occ_func_1_0(54)) + (occ_func_1_0(59)*occ_func_1_0(45)) + (occ_func_1_0(35)*occ_func_1_0(44)) + (occ_func_1_0(73)*occ_func_1_0(43)) + (occ_func_1_0(64)*occ_func_1_0(37)) + (occ_func_1_0(46)*occ_func_1_0(52)) + (occ_func_1_0(58)*occ_func_1_0(36)) + (occ_func_1_0(79)*occ_func_1_0(47)) + (occ_func_1_0(53)*occ_func_1_0(42)) + (occ_func_1_0(62)*occ_func_1_0(41)) + (occ_func_1_0(75)*occ_func_1_0(50)) + (occ_func_1_0(48)*occ_func_1_0(39)) + (occ_func_1_0(69)*occ_func_1_0(42)) + (occ_func_1_0(68)*occ_func_1_0(39)) + (occ_func_1_0(47)*occ_func_1_0(50)) + (occ_func_1_0(80)*occ_func_1_0(53)) + (occ_func_1_0(57)*occ_func_1_0(41)) + (occ_func_1_0(36)*occ_func_1_0(48)) + (occ_func_1_0(67)*occ_func_1_0(48)) + (occ_func_1_0(70)*occ_func_1_0(51)) + (occ_func_1_0(41)*occ_func_1_0(38)) + (occ_func_1_0(77)*occ_func_1_0(56)) + (occ_func_1_0(60)*occ_func_1_0(48)) + (occ_func_1_0(33)*occ_func_1_0(41)) + (occ_func_1_0(59)*occ_func_1_0(33)) + (occ_func_1_0(78)*occ_func_1_0(42)) + (occ_func_1_0(56)*occ_func_1_0(47)) + (occ_func_1_0(76)*occ_func_1_0(47)) + (occ_func_1_0(61)*occ_func_1_0(38)) + (occ_func_1_0(42)*occ_func_1_0(51)) + (occ_func_1_0(66)*occ_func_1_0(44)) + (occ_func_1_0(71)*occ_func_1_0(49)) + (occ_func_1_0(45)*occ_func_1_0(40)) + (occ_func_1_0(79)*occ_func_1_0(55)) + (occ_func_1_0(58)*occ_func_1_0(44)) + (occ_func_1_0(34)*occ_func_1_0(45)) + (occ_func_1_0(70)*occ_func_1_0(37)) + (occ_func_1_0(67)*occ_func_1_0(34)) + (occ_func_1_0(52)*occ_func_1_0(55)) + (occ_func_1_0(72)*occ_func_1_0(46)) + (occ_func_1_0(65)*occ_func_1_0(40)) + (occ_func_1_0(43)*occ_func_1_0(49)) + (occ_func_1_0(60)*occ_func_1_0(34)) + (occ_func_1_0(77)*occ_func_1_0(43)) + (occ_func_1_0(55)*occ_func_1_0(46)) + (occ_func_1_0(75)*occ_func_1_0(40)) + (occ_func_1_0(62)*occ_func_1_0(35)) + (occ_func_1_0(49)*occ_func_1_0(54)) + (occ_func_1_0(74)*occ_func_1_0(39)) + (occ_func_1_0(63)*occ_func_1_0(33)) + (occ_func_1_0(50)*occ_func_1_0(56)))/48.0;
  }

  /**** Basis functions for orbit 3, 20****
#Points: 3
MaxLength: 9.6513440  MinLength: 9.6513440
               0.5000000    0.5000000    0.5000000 Mn Ni
              -0.5000000   -0.5000000   -1.5000000 Mn Ni
               1.5000000   -1.5000000   -0.5000000 Mn Ni
****/
  double La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_20_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(69)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(79)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(73)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(57)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(70)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(75)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(61)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(74)))/8.0;
  }

  double La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_20_0() const{
    return ((occ_func_1_0(0)*occ_func_1_0(63)*occ_func_1_0(69)) + (occ_func_1_0(74)*occ_func_1_0(0)*occ_func_1_0(78)) + (occ_func_1_0(68)*occ_func_1_0(59)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(79)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(76)) + (occ_func_1_0(58)*occ_func_1_0(61)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(67)*occ_func_1_0(73)) + (occ_func_1_0(70)*occ_func_1_0(0)*occ_func_1_0(77)) + (occ_func_1_0(64)*occ_func_1_0(60)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(57)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(62)) + (occ_func_1_0(80)*occ_func_1_0(75)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(64)*occ_func_1_0(70)) + (occ_func_1_0(73)*occ_func_1_0(0)*occ_func_1_0(77)) + (occ_func_1_0(67)*occ_func_1_0(60)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(65)*occ_func_1_0(75)) + (occ_func_1_0(72)*occ_func_1_0(0)*occ_func_1_0(80)) + (occ_func_1_0(62)*occ_func_1_0(57)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(71)*occ_func_1_0(61)) + (occ_func_1_0(66)*occ_func_1_0(0)*occ_func_1_0(58)) + (occ_func_1_0(76)*occ_func_1_0(79)*occ_func_1_0(0)) + (occ_func_1_0(0)*occ_func_1_0(68)*occ_func_1_0(74)) + (occ_func_1_0(69)*occ_func_1_0(0)*occ_func_1_0(78)) + (occ_func_1_0(63)*occ_func_1_0(59)*occ_func_1_0(0)))/8.0;
  }

  double La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_20_0(int occ_i, int occ_f) const{
    return (m_occ_func_1_0[occ_f] - m_occ_func_1_0[occ_i])*((occ_func_1_0(63)*occ_func_1_0(69)) + (occ_func_1_0(74)*occ_func_1_0(78)) + (occ_func_1_0(68)*occ_func_1_0(59)) + (occ_func_1_0(71)*occ_func_1_0(79)) + (occ_func_1_0(66)*occ_func_1_0(76)) + (occ_func_1_0(58)*occ_func_1_0(61)) + (occ_func_1_0(67)*occ_func_1_0(73)) + (occ_func_1_0(70)*occ_func_1_0(77)) + (occ_func_1_0(64)*occ_func_1_0(60)) + (occ_func_1_0(65)*occ_func_1_0(57)) + (occ_func_1_0(72)*occ_func_1_0(62)) + (occ_func_1_0(80)*occ_func_1_0(75)) + (occ_func_1_0(64)*occ_func_1_0(70)) + (occ_func_1_0(73)*occ_func_1_0(77)) + (occ_func_1_0(67)*occ_func_1_0(60)) + (occ_func_1_0(65)*occ_func_1_0(75)) + (occ_func_1_0(72)*occ_func_1_0(80)) + (occ_func_1_0(62)*occ_func_1_0(57)) + (occ_func_1_0(71)*occ_func_1_0(61)) + (occ_func_1_0(66)*occ_func_1_0(58)) + (occ_func_1_0(76)*occ_func_1_0(79)) + (occ_func_1_0(68)*occ_func_1_0(74)) + (occ_func_1_0(69)*occ_func_1_0(78)) + (occ_func_1_0(63)*occ_func_1_0(59)))/8.0;
  }

}


extern "C" {
  /// \brief Returns a Clexulator_impl::Base* owning a La_MnNi_O3_cubic_Clexulator
  CASM::Clexulator_impl::Base* make_La_MnNi_O3_cubic_Clexulator() {
    return new CASM::La_MnNi_O3_cubic_Clexulator();
  }

}

