# Copyright 2019-2018 Jan-Niclas Luy
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import sys
import numpy as np

from ase.io import read as ase_read

from scipy.constants import physical_constants as pc

from nomadcore.simple_parser import SimpleMatcher
from nomadcore.baseclasses import ParserInterface, MainHierarchicalParser

from nomad.parsing import LocalBackend


"""
A very basic BAND parser.
"""


class BANDParser(ParserInterface):

    def get_metainfo_filename(self):
        return 'band.nomadmetainfo.json'

    def get_parser_info(self):
        return {
            'name': 'band_parser',
            'version': '1.0.0'
        }

    def setup_version(self):
        self.setup_main_parser(None)

    def setup_main_parser(self, _):
        self.main_parser = MainParser(self.parser_context)


class MainParser(MainHierarchicalParser):
    def __init__(self, parser_context, *args, **kwargs):
        super().__init__(parser_context, *args, **kwargs)
        self.lattice_vectors = []
        self.atom_labels = []
        self.atom_positions = []
        self.configuration_periodic_dimensions = []
        self.XC_functional_name = []
        self.LDA_functional_name = []
        self.GGA_functional_name = []
        self.meta_GGA_functional_name = []
        self.dos_values = []
        self.dos_energies = []

        self.root_matcher = SimpleMatcher(
            name='root',
            #adHoc=self.hallo,
            startReStr=r'B A N D',
            weak=True,
            sections=['section_run', 'section_method', 'section_system'],
            subMatchers=[
                SimpleMatcher(
                    startReStr=r' *   Amsterdam Density Functional  \(ADF\)\s*2018\s*\*.*',
                    subMatchers=[
                        SimpleMatcher(startReStr=r'\s\*\s{47}r(?P<program_version>\d+).*')
                    ],    
                    endReStr=r'\s\*{2}.*'),
                SimpleMatcher(
                    startReStr=r'Geometry.*',# G E O M E T R Y    I N    X - Y - Z    F O R M A T
                    subMatchers=[
                        SimpleMatcher(         
                            startReStr=r'  Index Symbol       x \(bohr\)       y \(bohr\)       z \(bohr\).*',
                            subMatchers=[
                                SimpleMatcher(startReStr=r'\s*\d+\s+([A-Z][a-z]?)\s+(\d+\.\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)\s*', repeats=True, startReAction=self.save_atoms)#might need to accoutn for negative sign
                            ],
                            endReStr=r'\s*'
                        ),
                        SimpleMatcher(
                            startReStr=r'Lattice vectors \(bohr\).*',#   REAL SPACE LATTICE VECTORS
                            subMatchers=[
                                SimpleMatcher(startReStr=r'\s*\d+\s+(\d+\.\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)\s*', repeats=True, startReAction=self.save_lattice)
                            ],
                            endReStr=r'\s*'
                        )                            
                    ]    
                ),
                SimpleMatcher(startReStr=r'Total System Charge             (?P<total_charge>\d+)\.\d+.*',
                ),
                SimpleMatcher( 
                     startReStr=r'Band Engine Input.*',
                     subMatchers=[
                         SimpleMatcher(
                             startReStr=r'((?i) \s*Basis)',
                             subMatchers=[
                                 SimpleMatcher(startReStr=r'((?i)\s*Type\s*(?P<basis_set>[a-z,A-Z]*).*)', )
                             ],
                             endReStr=r'((?i) \s*End)'
                         )
                    ],
                    endReStr=r' Using the following basis set files:'
                ),                               
                SimpleMatcher(
                    #sections=['section_XC_functionals'],
                    startReStr=r' DENSITY FUNCTIONAL POTENTIAL \(scf\)',
                    subMatchers=[
                        SimpleMatcher(startReStr=r'    LDA:                               ([a-z,A-Z]*).*', startReAction=self.save_lda
                        ),
                        SimpleMatcher(startReStr=r'    Gradient Corrections:              ([a-z,A-Z]*)c\s*([a-z,A-Z]*)x.*', startReAction=self.save_functional
                        ),
                        SimpleMatcher(startReStr=r'    Meta-GGA:                          ([a-z,A-Z]*).*', startReAction=self.save_meta_gga
                        )                        
                    ],
                    endReStr=r' DENSITY FUNCTIONAL ENERGY \(post-scf\)'
                ),
                SimpleMatcher(
                    startReStr=r' R U N    C O N F I G.*',
                    sections=['section_single_configuration_calculation'],
                    subMatchers=[
                        SimpleMatcher(
                            sections=['section_dos'],
                            startReStr=r' TOTALDOS.*',
                            subMatchers=[
                                SimpleMatcher(
                                    startReStr=r'\s*([+-]?\d+\.\d+E[+-]?\d+)\s*([+-]?\d+\.\d+E[+-]?\d+)\s*([+-]?\d+\.\d+E[+-]?\d+).*', repeats=True, startReAction=self.save_dos)
                            ],
                            endReStr=r' ENDINPUT.*',
                        ),
                        SimpleMatcher(startReStr=r'Energy \(hartree\)            (?P<energy_total__hartree>-\d+\.\d+).*'
                        )   
                    ]
                ),                             
            ]
        )
        
    def hallo(*args):
        print("hallo")
        
    def save_dos(self, _, groups):
        self.dos_values.append([float(groups[1]), float(groups[2])])
        #print('test')
        self.dos_energies.append(float(groups[0]))
             
    def save_atoms(self, _, groups):
        self.atom_positions.append([float(groups[1]), float(groups[2]), float(groups[3])])
        self.atom_labels.append(groups[0])
        
    def save_lattice(self, _, groups):
        self.lattice_vectors.append([float(groups[0]), float(groups[1]), float(groups[2])])

    def save_lda(self, _, groups):
        self.LDA_functional_name.append(groups[0])
        
    def save_meta_gga(self, _, groups):
        self.meta_GGA_functional_name.append(groups[0])
            
    def save_functional(self, _, groups):
        if groups != None:
            self.GGA_functional_name.append([groups[0], groups[1]])
        
    def onClose_section_system(self, backend, *args, **kwargs):
        backend.addArrayValues('atom_labels', np.array(self.atom_labels))
        backend.addArrayValues('atom_positions', np.array(self.atom_positions)*pc['Bohr radius'][0])
        #print(self.lattice_vectors)
        for _ in range(0, len(self.lattice_vectors)):
            self.configuration_periodic_dimensions.append(True)
        for _ in range(len(self.lattice_vectors),3):
            self.configuration_periodic_dimensions.append(False)
        #print(self.configuration_periodic_dimensions)
        for _ in range(len(self.lattice_vectors),3):
            self.lattice_vectors.append([0,0,0])
        backend.addArrayValues('lattice_vectors', np.array(self.lattice_vectors)*pc['Bohr radius'][0])
        backend.addArrayValues('configuration_periodic_dimensions', np.array(self.configuration_periodic_dimensions))
        
    def onClose_section_dos(self, backend, *args, **kwargs):
        #print(self.dos_energies)
        backend.addArrayValues('dos_values', np.array(self.dos_values))
        backend.addArrayValues('dos_energies', np.array(self.dos_energies)*pc['Hartree energy'][0])

    def onClose_section_method(self, backend, *args, **kwargs):
        backend.addValue('electronic_structure_method', 'DFT')
        if self.GGA_functional_name != []:
            #print(self.GGA_functional_name)
            backend.openNonOverlappingSection('section_XC_functionals')
            backend.addValue('XC_functional_name', 'GGA_X_' + self.GGA_functional_name[0][1])
            backend.closeNonOverlappingSection('section_XC_functionals')
            backend.openNonOverlappingSection('section_XC_functionals')
            backend.addValue('XC_functional_name', 'GGA_C_' + self.GGA_functional_name[0][0])
            backend.closeNonOverlappingSection('section_XC_functionals')
        if self.meta_GGA_functional_name != []:
            #print(self.meta_GGA_functional_name)
            backend.openNonOverlappingSection('section_XC_functionals')
            backend.addValue('XC_functional_name', 'MGGA_XC_' + self.meta_GGA_functional_name[0])
            backend.closeNonOverlappingSection('section_XC_functionals')
        if self.LDA_functional_name != [] and self.GGA_functional_name == []:
            #print(self.LDA_functional_name)
            backend.openNonOverlappingSection('section_XC_functionals')
            backend.addValue('XC_functional_name', 'LDA_XC_' + self.LDA_functional_name[0])
            backend.closeNonOverlappingSection('section_XC_functionals')

    def onClose_section_run(self, backend, *args, **kwargs):
        backend.addValue('program_name', 'BAND')
        backend.addValue('program_basis_set_type', 'numeric AOs')

if __name__ == "__main__":
    parser = BANDParser(backend=LocalBackend)
    parser.parse(sys.argv[1])
    #parser.parser_context.super_backend.write_json(sys.stdout)
    print(parser.parser_context.super_backend)
    
#TODO
# .out?
#geoopt?
