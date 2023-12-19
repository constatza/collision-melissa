import logging
import random

from melissa.server.sensitivity_analysis import SensitivityAnalysisServer
from melissa.server.parameters import RandomUniform  # , HaltonGenerator
from typing import Dict, Any

logger = logging.getLogger("melissa")
random.seed(123)


class HeatPDEServerSA(SensitivityAnalysisServer):
    """
    Use-case specific server
    """

    def __init__(self, config: Dict[str, Any]):
        super().__init__(config)
        self.nb_parms = self.study_options['nb_parameters']
        Tmin, Tmax = self.study_options['parameter_range']
        nb_sim = self.study_options['parameter_sweep_size']
        # Example of random uniform sampling
        self.parameter_generator = RandomUniform(
            nb_parms=self.nb_parms, nb_sim=nb_sim,
            l_bounds=[Tmin], u_bounds=[Tmax],
            second_order=False,
            apply_pick_freeze=self.sobol_op
        ).generator()

        # Example of Halton sampling
        # self.parameter_generator = HaltonGenerator(nb_parms=self.nb_parms, nb_sim=nb_sim,
        #                                            l_bounds=[Tmin],u_bounds=[Tmax],second_order=False,
        #                                            apply_pick_freeze=self.sobol_op).generator()

        # Example of Latin Hypercube Sampling
        # self.parameter_generator = LHSGenerator(nb_parms=self.nb_parms, nb_sim=nb_sim,
        #                                            l_bounds=[Tmin],u_bounds=[Tmax],second_order=False,
        #                                            apply_pick_freeze=self.sobol_op).generator()
