from modules import refactor,methylation_data
from numpy import loadtxt, array_equal
from utils import LinearRegression
import logging
from tests import tools
from pickle import load

class FeatureSelectionTester():
    FAKE_DATA  = "tests/refactor/files/feature_selection/data"
    FAKE_CONTROL = "tests/refactor/files/feature_selection/control"
    def __init__(self):
        self.fs_meth_data = methylation_data.MethylationData(datafile = self.FAKE_DATA)
        self.test_controls_fs()
        self.test_phenotype_fs()

    def test_controls_fs(self):
        """
        test that all samples that are '1' in the FAKE_CONTROL file are deleted from FAKE_DATA file

        """
        logging.info("Testing controls feature selection...")
        refactor_meth_data = self.fs_meth_data.copy()

        controls = loadtxt(self.FAKE_CONTROL, dtype = str)
        controls_indices = controls[:,1]
        module  = refactor.Refactor(methylation_data = refactor_meth_data, 
                      k = 2, 
                      feature_selection = "controls",
                      phenofile = self.FAKE_CONTROL,
                      t=5)
        index =0
        res = module.feature_selection_handler()

        # controls fs suppose to change methylation data
        assert array_equal(module.meth_data.data, res)

        # res is the data matrix found in FAKE_DATA without the columns (samples) at index i when i is index of a 0 in FAKE_CONTROL file
        for i in range(len(controls_indices)):
            if controls_indices[i]=='0':
                assert array_equal(self.fs_meth_data.data[:,i] ,res[:,index])
                index += 1

        # controls fs result will contain  number of columns (samples) as the number of 0's in the FAKE_CONTROL file
        assert index == len(res[0])

        
        logging.info("PASS")

    def test_phenotype_fs(self):
        logging.info("Testing phenotype feature selection...")
        refactor_meth_data = self.fs_meth_data.copy()

        module  = refactor.Refactor(methylation_data = refactor_meth_data, 
                      k = 2, 
                      feature_selection = "phenotype",
                      phenofile = self.FAKE_CONTROL,
                      t=5)

        phenotype = module._validate_phenotype(self.FAKE_CONTROL, "phenotype")
        # validate phenotype feature selection output (res_data) is correlated to our linear regression for (site, phenotype)
        res_data = module.feature_selection_handler()
        for i,site in enumerate(self.fs_meth_data.data):
            lin_reg = LinearRegression(site, phenotype)
            assert len(lin_reg.residuals) == len(site)
            # validate our residuals are corelated to res_data
            assert tools.correlation(lin_reg.residuals, res_data[i])

        logging.info("PASS")

class RefactorTester():
    DEMO_SMALL_DATA = "tests/refactor/files/demofiles/datafile2"
    DEMO_DATA_NO_BAD_PROBES = "tests/refactor/files/demofiles/datafile2_no_bad_probes"
    DEMO_COVAR = "tests/refactor/files/demofiles/covariates"
    DEMO_PHENO = "tests/refactor/files/demofiles/phenotype"
    DEMO_CELLPRO = "tests/refactor/files/demofiles/cellproportions"
    BAD_PROBES = "tests/refactor/files/demofiles/bad_probes"

    # senarios output
    COMP_K5_T400 = "tests/refactor/files/senarios_out/k5t400.out.components.txt"
    RANK_K5_T400 = "tests/refactor/files/senarios_out/k5t400.out.rankedlist.txt"
    
    COMP_K5_T400_stdth01numcomp7 = "tests/refactor/files/senarios_out/k5t400stdth0.1numcomp7.out.components.txt"
    RANK_K5_T400_stdth01numcomp7 = "tests/refactor/files/senarios_out/k5t400stdth0.1numcomp7.out.rankedlist.txt"
    
    COMP_K5_T400_stdth013 = "tests/refactor/files/senarios_out/k5t400stdth0.13.out.components.txt"
    RANK_K5_T400_stdth013 = "tests/refactor/files/senarios_out/k5t400stdth0.13.out.rankedlist.txt"
    
    COMP_K5_T400_covar = "tests/refactor/files/senarios_out/k5t400covar.out.components.txt"
    RANK_K5_T400_covar = "tests/refactor/files/senarios_out/k5t400covar.out.rankedlist.txt"
  
    COMP_K5_T400_stdth008covar = "tests/refactor/files/senarios_out/k5t400stdth0.08covar.out.components.txt"
    RANK_K5_T400_stdth008covar = "tests/refactor/files/senarios_out/k5t400stdth0.08covar.out.rankedlist.txt"

    def __init__(self):
        self.meth_data = methylation_data.MethylationData(datafile = self.DEMO_SMALL_DATA)
        self.test_remove_covariates()
        self.test_low_rank_approx_distances()
        self.test_exclude_bad_probes()
        self.test_senarios()

    def test_remove_covariates(self):
        logging.info("Testing removing covariates...")
        covar_meth_data = self.meth_data.copy()

        module  = refactor.Refactor(methylation_data = covar_meth_data, 
                      k = 2, 
                      covar = self.DEMO_COVAR,
                      t=500)

        coavr = module._validate_covar(self.DEMO_COVAR,)
        # remove from refactor_meth_data
        module._remove_covariates()

        # remove "manually"
        for i,site in enumerate(self.meth_data.data):
            lin_reg = LinearRegression(site, coavr)
            assert len(lin_reg.residuals) == len(site)

            # validate our residuals are corelated to res_data
            assert tools.correlation(lin_reg.residuals, covar_meth_data.data[i])

        logging.info("PASS")


    def test_low_rank_approx_distances(self):
        """
        tests that number of distances is as the number of sites (distance for every site)
        """
        logging.info("Testing low rank approx distances...")
        dis_meth_data = self.meth_data.copy()

        module  = refactor.Refactor(methylation_data = dis_meth_data, 
                                    k = 5)

        distances = module._calc_low_rank_approx_distances()
        assert distances.size == dis_meth_data.sites_size, "there must be distances as the number of sites"
        logging.info("PASS")

    
    def test_exclude_bad_probes(self):
        logging.info("Testing removing bad probes...")
        probes_meth_data = self.meth_data.copy()

        data_no_bad_probes = methylation_data.MethylationData(datafile = self.DEMO_DATA_NO_BAD_PROBES)

        bad_probes = load(open(self.BAD_PROBES,'r'))
        module  = refactor.Refactor(methylation_data = probes_meth_data, 
                                    k = 5, 
                                    bad_probes_list = bad_probes)

        module._exclude_bad_probes()


        assert array_equal(data_no_bad_probes.data, module.meth_data.data)

        # tests sites list has changed
        remove_count = len(bad_probes)
        orig_sites_before = []
        orig_sites_before.extend(self.meth_data.cpgnames)
        orig_sites_after = []
        orig_sites_after.extend(module.meth_data.cpgnames)
        for i in bad_probes:
            try:
                orig_sites_before.remove(i)
            except:
                remove_count -= 1
        assert orig_sites_after == orig_sites_before
        # test sites size
        assert self.meth_data.sites_size - remove_count == module.meth_data.sites_size

        logging.info("PASS")

    def test_senarios(self):

        logging.info("Testing senario no.1...")
        senario_meth_data = self.meth_data.copy()
        module  = refactor.Refactor(methylation_data = senario_meth_data, 
                                    k = 5, 
                                    t = 400)
        module.run()
        comp = loadtxt(self.COMP_K5_T400)
        assert module.components.shape == comp.shape

        for i in range(module.components.shape[1]):
            assert tools.correlation(module.components[:,i], comp[:,i])

        assert array_equal(module.ranked_sites, loadtxt(self.RANK_K5_T400, dtype = str))
        logging.info("PASS")

        logging.info("Testing senario no.2...")
        senario_meth_data = self.meth_data.copy()
        module  = refactor.Refactor(methylation_data = senario_meth_data, 
                                      k = 5, 
                                      t = 400,
                                      minstd = 0.1,
                                      num_components = 7)
        module.run()
        comp = loadtxt(self.COMP_K5_T400_stdth01numcomp7)
        assert module.components.shape == comp.shape

        for i in range(module.components.shape[1]):
            assert tools.correlation(module.components[:,i], comp[:,i])

        assert array_equal(module.ranked_sites, loadtxt(self.RANK_K5_T400_stdth01numcomp7, dtype = str))
        logging.info("PASS")

        logging.info("Testing senario no.3...")
        senario_meth_data = self.meth_data.copy()
        module  = refactor.Refactor(methylation_data = senario_meth_data, 
                                      k = 5, 
                                      t = 400,
                                      minstd = 0.13)
        module.run()
        comp = loadtxt(self.COMP_K5_T400_stdth013)
        assert module.components.shape == comp.shape

        for i in range(module.components.shape[1]):
            assert tools.correlation(module.components[:,i], comp[:,i])

        assert array_equal(module.ranked_sites, loadtxt(self.RANK_K5_T400_stdth013, dtype = str))
        logging.info("PASS")



        logging.info("Testing senario no.4...")
        senario_meth_data = self.meth_data.copy()
        module  = refactor.Refactor(methylation_data = senario_meth_data, 
                                      k = 5, 
                                      t = 400,
                                      covar = self.DEMO_COVAR)
        module.run()
        comp = loadtxt(self.COMP_K5_T400_covar)
        assert module.components.shape == comp.shape

        for i in range(module.components.shape[1]):
            assert tools.correlation(module.components[:,i], comp[:,i])

        assert array_equal(module.ranked_sites, loadtxt(self.RANK_K5_T400_covar, dtype = str))
        logging.info("PASS")


        logging.info("Testing senario no.5...")
        senario_meth_data = self.meth_data.copy()
        module  = refactor.Refactor(methylation_data = senario_meth_data, 
                                      k = 5, 
                                      t = 400,
                                      minstd = 0.08,
                                      covar = self.DEMO_COVAR)
        module.run()

        comp = loadtxt(self.COMP_K5_T400_stdth008covar)
        assert module.components.shape == comp.shape

        for i in range(module.components.shape[1]):
            assert tools.correlation(module.components[:,i], comp[:,i])

        assert array_equal(module.ranked_sites, loadtxt(self.RANK_K5_T400_stdth008covar, dtype = str))
        logging.info("PASS")

  

