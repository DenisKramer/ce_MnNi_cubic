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
    BasisFuncPtr m_orbit_func_list[10];

    // array of pointers to member functions for calculating flower functions
    BasisFuncPtr m_flower_func_lists[5][10];

    // array of pointers to member functions for calculating DELTA flower functions
    DeltaBasisFuncPtr m_delta_func_lists[5][10];

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

    double eval_bfunc_3_0_0() const;

    double site_eval_at_1_bfunc_3_0_0() const;

    double delta_site_eval_at_1_bfunc_3_0_0(int occ_i, int occ_f) const;

    double eval_bfunc_3_1_0() const;

    double site_eval_at_1_bfunc_3_1_0() const;

    double delta_site_eval_at_1_bfunc_3_1_0(int occ_i, int occ_f) const;


  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  La_MnNi_O3_cubic_Clexulator::La_MnNi_O3_cubic_Clexulator() :
    Clexulator_impl::Base(81, 10) {
    m_occ_func_1_0[0] = -0.0000000000, m_occ_func_1_0[1] = 1.0000000000;

    m_orbit_func_list[0] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_0_0_0;
    m_orbit_func_list[1] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_1_0_0;
    m_orbit_func_list[2] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_0_0;
    m_orbit_func_list[3] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_1_0;
    m_orbit_func_list[4] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_2_0;
    m_orbit_func_list[5] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_3_0;
    m_orbit_func_list[6] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_4_0;
    m_orbit_func_list[7] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_2_5_0;
    m_orbit_func_list[8] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_0_0;
    m_orbit_func_list[9] = &La_MnNi_O3_cubic_Clexulator::eval_bfunc_3_1_0;


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


    m_flower_func_lists[1][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_flower_func_lists[1][1] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_1_0_0;
    m_flower_func_lists[1][2] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_0_0;
    m_flower_func_lists[1][3] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_1_0;
    m_flower_func_lists[1][4] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_2_0;
    m_flower_func_lists[1][5] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_3_0;
    m_flower_func_lists[1][6] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_4_0;
    m_flower_func_lists[1][7] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_2_5_0;
    m_flower_func_lists[1][8] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_0_0;
    m_flower_func_lists[1][9] = &La_MnNi_O3_cubic_Clexulator::site_eval_at_1_bfunc_3_1_0;


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


    m_delta_func_lists[1][0] = &La_MnNi_O3_cubic_Clexulator::zero_func;
    m_delta_func_lists[1][1] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_1_0_0;
    m_delta_func_lists[1][2] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_0_0;
    m_delta_func_lists[1][3] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_1_0;
    m_delta_func_lists[1][4] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_2_0;
    m_delta_func_lists[1][5] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_3_0;
    m_delta_func_lists[1][6] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_4_0;
    m_delta_func_lists[1][7] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_2_5_0;
    m_delta_func_lists[1][8] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_0_0;
    m_delta_func_lists[1][9] = &La_MnNi_O3_cubic_Clexulator::delta_site_eval_at_1_bfunc_3_1_0;


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


    m_weight_matrix.row(0) << 1, 0, 0;
    m_weight_matrix.row(1) << 0, 1, 0;
    m_weight_matrix.row(2) << 0, 0, 1;

    m_neighborhood = std::set<UnitCellCoord> {
      {UnitCellCoord(1, -2, -1, -1)},
      {UnitCellCoord(1, -2, -1, 0)},
      {UnitCellCoord(1, -2, -1, 1)},
      {UnitCellCoord(1, -2, 0, -1)},
      {UnitCellCoord(1, -2, 0, 0)},
      {UnitCellCoord(1, -2, 0, 1)},
      {UnitCellCoord(1, -2, 1, -1)},
      {UnitCellCoord(1, -2, 1, 0)},
      {UnitCellCoord(1, -2, 1, 1)},
      {UnitCellCoord(1, -1, -2, -1)},
      {UnitCellCoord(1, -1, -2, 0)},
      {UnitCellCoord(1, -1, -2, 1)},
      {UnitCellCoord(1, -1, -1, -2)},
      {UnitCellCoord(1, -1, -1, -1)},
      {UnitCellCoord(1, -1, -1, 0)},
      {UnitCellCoord(1, -1, -1, 1)},
      {UnitCellCoord(1, -1, -1, 2)},
      {UnitCellCoord(1, -1, 0, -2)},
      {UnitCellCoord(1, -1, 0, -1)},
      {UnitCellCoord(1, -1, 0, 0)},
      {UnitCellCoord(1, -1, 0, 1)},
      {UnitCellCoord(1, -1, 0, 2)},
      {UnitCellCoord(1, -1, 1, -2)},
      {UnitCellCoord(1, -1, 1, -1)},
      {UnitCellCoord(1, -1, 1, 0)},
      {UnitCellCoord(1, -1, 1, 1)},
      {UnitCellCoord(1, -1, 1, 2)},
      {UnitCellCoord(1, -1, 2, -1)},
      {UnitCellCoord(1, -1, 2, 0)},
      {UnitCellCoord(1, -1, 2, 1)},
      {UnitCellCoord(1, 0, -2, -1)},
      {UnitCellCoord(1, 0, -2, 0)},
      {UnitCellCoord(1, 0, -2, 1)},
      {UnitCellCoord(1, 0, -1, -2)},
      {UnitCellCoord(1, 0, -1, -1)},
      {UnitCellCoord(1, 0, -1, 0)},
      {UnitCellCoord(1, 0, -1, 1)},
      {UnitCellCoord(1, 0, -1, 2)},
      {UnitCellCoord(1, 0, 0, -2)},
      {UnitCellCoord(1, 0, 0, -1)},
      {UnitCellCoord(1, 0, 0, 0)},
      {UnitCellCoord(1, 0, 0, 1)},
      {UnitCellCoord(1, 0, 0, 2)},
      {UnitCellCoord(1, 0, 1, -2)},
      {UnitCellCoord(1, 0, 1, -1)},
      {UnitCellCoord(1, 0, 1, 0)},
      {UnitCellCoord(1, 0, 1, 1)},
      {UnitCellCoord(1, 0, 1, 2)},
      {UnitCellCoord(1, 0, 2, -1)},
      {UnitCellCoord(1, 0, 2, 0)},
      {UnitCellCoord(1, 0, 2, 1)},
      {UnitCellCoord(1, 1, -2, -1)},
      {UnitCellCoord(1, 1, -2, 0)},
      {UnitCellCoord(1, 1, -2, 1)},
      {UnitCellCoord(1, 1, -1, -2)},
      {UnitCellCoord(1, 1, -1, -1)},
      {UnitCellCoord(1, 1, -1, 0)},
      {UnitCellCoord(1, 1, -1, 1)},
      {UnitCellCoord(1, 1, -1, 2)},
      {UnitCellCoord(1, 1, 0, -2)},
      {UnitCellCoord(1, 1, 0, -1)},
      {UnitCellCoord(1, 1, 0, 0)},
      {UnitCellCoord(1, 1, 0, 1)},
      {UnitCellCoord(1, 1, 0, 2)},
      {UnitCellCoord(1, 1, 1, -2)},
      {UnitCellCoord(1, 1, 1, -1)},
      {UnitCellCoord(1, 1, 1, 0)},
      {UnitCellCoord(1, 1, 1, 1)},
      {UnitCellCoord(1, 1, 1, 2)},
      {UnitCellCoord(1, 1, 2, -1)},
      {UnitCellCoord(1, 1, 2, 0)},
      {UnitCellCoord(1, 1, 2, 1)},
      {UnitCellCoord(1, 2, -1, -1)},
      {UnitCellCoord(1, 2, -1, 0)},
      {UnitCellCoord(1, 2, -1, 1)},
      {UnitCellCoord(1, 2, 0, -1)},
      {UnitCellCoord(1, 2, 0, 0)},
      {UnitCellCoord(1, 2, 0, 1)},
      {UnitCellCoord(1, 2, 1, -1)},
      {UnitCellCoord(1, 2, 1, 0)},
      {UnitCellCoord(1, 2, 1, 1)}
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

    m_orbit_neighborhood[9] = std::set<UnitCellCoord> {
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

}


extern "C" {
  /// \brief Returns a Clexulator_impl::Base* owning a La_MnNi_O3_cubic_Clexulator
  CASM::Clexulator_impl::Base* make_La_MnNi_O3_cubic_Clexulator() {
    return new CASM::La_MnNi_O3_cubic_Clexulator();
  }

}

