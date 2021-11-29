from get_testcase_data import get_testcase_data
from create_ic import create_ic
from run_model import run_model
from stress_divergence_map import stress_divergence_map
from stress_divergence_scaling import stress_divergence_scaling
from error_analysis_stress_divergence import error_analysis_stress_divergence

#-------------------------------------------------------------------------------

def run_stress_divergence_testcase():

    get_testcase_data()

    create_ic()

    run_model()

    #stress_divergence_map()

    stress_divergence_scaling()
     
    error_analysis_stress_divergence()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_stress_divergence_testcase()
