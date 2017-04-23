#!/usr/bin/env python
import sys, os, pandas, numpy, pylab
import traceback
import sbml_mca


class SBMLTestSuite:
    """
    class for running the SBML test suite (sbml.org/Software/SBML_Test_Suite)

    date:   04. July 2014
    author: Jannis Uhlendorf
    """

    _accepted_errors = ['Error: Events not supported yet',\
                       'Empty stoichiometric matrix.',\
                       'Initial assignments are currently only supported for parameters.',\
                       'Error: Varying compartment sizes not yet supported',\
                       'Algebraic rules not supported' ]
    _TRESHOLD = 5e-2

    def __init__(self, path, writer=sys.stdout):
        """
        class constructor
        @param path: path to the SBML test suite (the directory containing the cases directory)
        @param writer: (optional) object implementing writer (default  sys.stdout)
        """
        self._path    =  os.path.join( path + '/cases/semantic' )
        print self._path
        self._writer  =  writer

    def run_test( self, test_no = None, stop_and_plot=False ):
        """
        run either a single or all tests.
        if test_no is specified, only this test is done, if test_no is empy, all tests are done
        @param test_no: integer specifiyng which test to run
        """
        no_passed           = 0
        no_failed_distance  = 0
        no_failed_exception = 0
        no_accepted_errors  = 0
        no_unknown_errors   = 0
        unknown_errors      = {}
        # generate list of test case numbers
        if test_no == None:
            tests = range( 1, 1196 )
        elif isinstance(test_no, list):
            tests = test_no
        else:
            tests = [test_no]            

        for number in tests:
            self._print( 'Running test no. ' + str(number) )
            if number in [1108]: #[327,328,329,331,332,333,334]:
                # for these cases, something went seriously wrong, so for now we skip them
                no_failed_exception += 1
                continue

            # generate SBML path
            sbml_version_prioriy = ('l2v4', 'l2v3', 'l2v2', 'l2v1', 'l3v1', 'l1v2') 
            f_path = None
            for ver in sbml_version_prioriy:
                p = self._generate_file_path( number, '-sbml-%s.xml'%ver )
                if os.path.isfile(p):
                    f_path = p
                    break
            if not f_path:
                print f_path
                self._print( "Error: no suitable SBML file was found for %d" %number )
                sys.exit()
                
            try:
                mca                = sbml_mca.sbml_mca( f_path )
                settings           = self._read_settings( number )
                reference_result   = self._read_results( number )

                d = self._run_single_test( mca, settings, reference_result )
                if d > self._TRESHOLD * len(settings['variables']):
                    no_failed_distance += 1
                    self._print( '\t FAILED with simulation distance %e' %d )
                    if stop_and_plot:
                        self._plot_simulation( mca, settings, reference_result )
                        sys.exit()
                else:
                    no_passed += 1
                    self._print( '\t OK - passed with simulation distance %e' %d )
            except Exception, e:
                if e.message in self._accepted_errors:
                    no_accepted_errors += 1
                else:
                    no_unknown_errors += 1
                    if unknown_errors.has_key(e.message):
                        unknown_errors[e.message] += 1
                    else:
                        unknown_errors[e.message] = 1
                #print traceback.format_exc()
                no_failed_exception += 1
                self._print( '\t FAILED with following error:' )
                self._print( '\t' + str(type(e)) + ': ' + e.message )

        self._print( '\n\n*** Finished testing ***')
        self._print( '\n %d of %d tests passed' %(no_passed, len(tests)))
        self._print( '\n %d tests failed because of wrong simulation results' %(no_failed_distance) )
        self._print( '\n %d tests failed because of exception (not supported)' %(no_failed_exception) )
        self._print( '\n %d accepted exceptions (we know what went wrong)' %(no_accepted_errors) )
        self._print( '\n %d unknown exceptions (we don\'t know what went wrong)' %(no_unknown_errors) )
        if unknown_errors != {}:
            self._print( '\n Unkown Errors that have occurred:')
            for key in unknown_errors:
                self._print( '\n\t' + key + ':\t'+ str(unknown_errors[key]) )
        self._print( '\n' )
                     


    def _run_single_test( self, mca_obj, settings, reference_result ):
        """
        run a single test
        @param mca_obj: sbml_mca object of the model to compare
        @param settings: settigns dictionary as returned by _read_settings fct.
        @param reference_results: result dictionary as returned by _compare_results fct.
        """
        #variables = mca_obj._species_ids + [s.getId() for s in mca_obj._get_constant_species()]
        #time, t_courses = self._simulate( mca_obj, settings )
        simulation_result = self._simulate( mca_obj, settings )
        #simulation_result = { 'time': time }
        #for pos,v in enumerate(variables):
        #    print v
        #    simulation_result[v] = t_courses[:,pos]
        return self._compare_results( reference_result, simulation_result )

    def _simulate( self, mca_obj, settings ):
        # simlate model with given settings
        #time, t_courses = mca_obj.integrate_with_constant_species( settings['duration'],
        #                                                           settings['steps']+1,
        #                                                           r_tol=settings['relative'],
        #                                                           a_tol=settings['absolute'] )
        result = mca_obj.integrate_return_dict( settings['duration'],
                                                settings['steps']+1,
                                                r_tol=settings['relative'],
                                                a_tol=settings['absolute'],
                                                with_constant_species=True,
                                                with_assignment_rules=True )
        return result

    def _print( self, msg ):
        """ print a message using the classes writer object """
        self._writer.write( str(msg) + '\n' )

    def _generate_file_path( self, number, ending ):
        number_string = '%05d' %number
        return os.path.join( self._path, number_string, number_string + ending )

    def _compare_results( self, result_a, result_b ):
        diff = 0.
        for key in result_a:
            if key=='time':
                continue
            d = numpy.linalg.norm( result_a[key] - result_b[key] )
            #print key, d
            diff += d
        return diff


    def _read_settings( self, number ):
        number_items = ['start', 'duration', 'steps', 'absolute', 'relative']
        list_items   = ['variables', 'amount', 'concentration']
        result = {}
        f_path = self._generate_file_path( number, '-settings.txt' )
        f = open( f_path, 'r' )
        for line in f.readlines():
            try:
                key, data = line.split(': ')
            except:
                pass
                #print line

            if key in number_items:
                result[key] = float(data)
            elif key in list_items:
                try:
                    result[key] = [ x.strip() for x  in data.split(', ' ) ]
                except:
                    result[key] = data.strip()
            else:
                result[key] = data
        f.close()
        return result

    def _read_results( self, number ):
        f_path = self._generate_file_path( number, '-results.csv' )
        result = {}
        csv_result = pandas.read_csv( f_path )
        variables = csv_result.axes[1]
        for v in variables:
            result[v] = numpy.array( csv_result[v].tolist() )
        return result

    def _plot_simulation( self, mca_obj, settings, reference_result ):
        variables = mca_obj._species_ids + [s.getId() for s in mca_obj._get_constant_species()]
        #time, t_courses = self._simulate( mca_obj, settings )
        simulation_result = self._simulate( mca_obj, settings )

        print reference_result
        print simulation_result
        
        colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
        for pos,v in enumerate(variables):
            pylab.plot( reference_result['time'], reference_result[v], ':x'+colors[pos] )
            pylab.plot( simulation_result['time'], simulation_result[v], '-'+colors[pos])
        pylab.show()

if __name__=='__main__':

    if len(sys.argv)<2:
        path = './sbml-test-cases'
    else:
        path = sys.argv[1]
    sts = SBMLTestSuite(path)
    #sts.run_test(  81, stop_and_plot=True)
    #sts.run_test(range(1,100))
    sts.run_test()
