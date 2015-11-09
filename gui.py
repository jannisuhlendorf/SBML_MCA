#!/usr/bin/env python


import enthought.traits.api
import enthought.traits.ui
import enthought.traits.ui.api
import enthought.traits.ui.menu
import enthought.traits.ui.file_dialog

#import enthought.chaco.chaco_plot_editor


import sbml_mca




class time_course(enthought.traits.api.HasTraits):
    sim_time = enthought.traits.api.Float(100, label="Simulation time", desc="End time of the simulation")
    steps    = enthought.traits.api.Float(100, label="Steps", desc="Number of time steps")
    run      = enthought.traits.api.Button('Run simulation')

    view = enthought.traits.ui.api.View( enthought.traits.ui.api.Item(name='sim_time'),
                                         enthought.traits.ui.api.Item(name='steps'),
                                         enthought.traits.ui.api.Item(name='run') )
    
    #def __init__( self,  **traits ):
        #self._super = traits['super']
        #del traits['super']
        #enthought.traits.api.HasTraits.__init__(self, **traits)

    def _set_sup( self, sup ):
        self._sup = sup
                  
    def _run_fired(self):
        self._sup.sim.plot_timecourse( self.sim_time, self.steps )


class Employee ( enthought.traits.api.HasTraits ):
    name  = enthought.traits.api.Str( '<unknown>' )
    title = enthought.traits.api.Str
    phone = enthought.traits.api.Regex( regex = r'\d\d\d-\d\d\d\d' )
class Department ( enthought.traits.api.HasTraits ):
    name      = enthought.traits.api.Str( '<unknown>' )
    employees = enthought.traits.api.List( Employee )
class Company ( enthought.traits.api.HasTraits ):
    name        = enthought.traits.api.Str( '<unknown>' )
    departments = enthought.traits.api.List( Department )
    employees   = enthought.traits.api.List( Employee )





class model(enthought.traits.api.HasTraits):


    class sbml_model(enthought.traits.api.HasTraits):
        name = enthought.traits.api.Str('')
        view = enthought.traits.ui.api.View( enthought.traits.ui.api.Item(name='name') )

    no_view = enthought.traits.ui.api.View()

    tree_editor = enthought.traits.ui.editors.TreeEditor( 
    nodes = [
        enthought.traits.ui.api.TreeNode( node_for  = [ Company ],
                  auto_open = True,
                  children  = '',
                  label     = 'name',
                  view      = enthought.traits.ui.api.View( [ 'name' ] )
        ),
        enthought.traits.ui.api.TreeNode( node_for  = [ Company ],
                  auto_open = True,
                  children  = 'departments',
                  label     = '=Departments',
                  view      = no_view,
                  add       = [ Department ]
        ),
        enthought.traits.ui.api.TreeNode( node_for  = [ Company ],
                  auto_open = True,
                  children  = 'employees',
                  label     = '=Employees',
                  view      = no_view,
                  add       = [ Employee ]
        ),
        enthought.traits.ui.api.TreeNode( node_for  = [ Department ],
                  auto_open = True,
                  children  = 'employees',
                  label     = 'name',
                  view      = enthought.traits.ui.api.View( [ 'name' ] ),
                  add       = [ Employee ]
        ),
        enthought.traits.ui.api.TreeNode( node_for  = [ Employee ],
                  auto_open = True,
                  label     = 'name',
                  view      = enthought.traits.ui.api.View( [ 'name', 'title', 'phone' ] )
                  )
        ] )

        
       
    def default_title ( self ):
        self.title = 'Senior Engineer'



    m = enthought.traits.api.Instance( sbml_model )


    model_tree  = enthought.traits.ui.editors.tree_editor.TreeEditor(nodes = [
        enthought.traits.ui.api.TreeNode( node_for  = [ sbml_model ],
                  auto_open = True,
                  children  = '',
                  label     = 'sdf',
                  view      = enthought.traits.ui.api.View(['name']) ) ] )

    view = enthought.traits.ui.api.View( enthought.traits.ui.api.Item(name="m", editor=model_tree, show_label=False ),
        title     = 'Company Structure',
        buttons   = [ 'OK' ],
        resizable = True,
        style     = 'custom',
        width     = .3,
        height    = .3 )
    

    
    


open_file = enthought.traits.ui.menu.Action(name = "Open File",  
                                            action = "open_file")

class main(enthought.traits.api.HasTraits):
    tc           =  enthought.traits.api.Instance( time_course, ())
    model_view   =  enthought.traits.api.Instance( model, () ) 
    
    def __init__(self, **traits):
        self.tc._set_sup(self)        
        enthought.traits.api.HasTraits.__init__(self, **traits)
        

    def load_file(self, f_name):
        self.sim = sbml_mca.sbml_mca( f_name )

        

    view = enthought.traits.ui.api.View( 
        enthought.traits.ui.api.Group( 
            enthought.traits.ui.api.Item("model_view", style='custom', show_label=False, dock="tab" ),
            enthought.traits.ui.api.Item("tc", style="custom", show_label=False, dock="tab"),
            layout="tabbed"),
        menubar = enthought.traits.ui.menu.MenuBar( enthought.traits.ui.menu.Menu( open_file, name='File' ) ),
        height=0.75, width=0.57 )



class main_handler(enthought.traits.ui.api.Controller):
    def open_file(self, info):
        f_name = str(enthought.traits.ui.file_dialog.open_file())
        if f_name:
            self.model.load_file( f_name )
        




#main_view = enthought.traits.ui.api.View( enthought.traits.ui.api.Item(name='first_name'),
#                                          enthought.traits.ui.api.Item(name='last_name'),
#                                          enthought.traits.ui.api.Item(name='department'),
#                                          enthought.traits.ui.api.Item(name='model'),
#                                          menubar = enthought.traits.ui.menu.MenuBar( enthought.traits.ui.menu.Menu( open_file, name='File' ) ),
#                                          height=0.75, width=0.57 )

gui = main()
gui.configure_traits( handler=main_handler(model=gui))
