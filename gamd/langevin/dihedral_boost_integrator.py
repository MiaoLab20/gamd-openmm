from abc import ABC

from gamd.langevin.base_integrator import GamdLangevinIntegrator


class GroupBoostIntegrator(GamdLangevinIntegrator, ABC):

    def __init__(self):
        pass

    def _add_common_variables(self):
        super(GamdLangevinIntegrator, self)._add_common_variables()

    def _add_conventional_md_pre_calc_step(self):
        pass

    def _add_conventional_md_position_update_step(self):
        pass

    def _add_conventional_md_velocity_update_step(self):
        pass

    def _add_conventional_md_stochastic_velocity_update_step(self):
        pass

    def _add_gamd_position_update_step(self):
        pass

    def _add_gamd_velocity_update_step(self):
        pass

    def _add_gamd_stochastic_velocity_update_step(self):
        pass

    def _add_gamd_pre_calc_step(self):
        pass

    def _add_gamd_boost_calculations_step(self):
        pass

    def _add_instructions_to_calculate_primary_boost_statistics(self):
        pass

    def _add_instructions_to_calculate_secondary_boost_statistics(self):
        pass

    def _update_potential_state_values_with_window_potential_state_values(self):
        pass

    def get_current_potential_energy(self):
        pass

    def get_total_force_scaling_factor(self):
        pass

    def get_boost_potential(self):
        pass


class DihedralBoostIntegrator(GroupBoostIntegrator, ABC):

    def __init__(self):
        pass


class LowerBoundIntegrator(DihedralBoostIntegrator):

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass


class UpperBoundIntegrator(DihedralBoostIntegrator):

    def _calculate_threshold_energy_and_effective_harmonic_constant(self):
        pass
