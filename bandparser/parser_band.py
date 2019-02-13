# Copyright 2019-2018 Markus Scheidgen
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
        self.GGA_functional_name = []

        self.root_matcher = SimpleMatcher(
            name='root',
			#adHoc=self.hallo,
            startReStr=r'B A N D',
            weak=True,
            sections=['section_run'],
            subMatchers=[
                SimpleMatcher(
				    startReStr=r'\s\*\s*Amsterdam\sDensity\sFunctional\s*\(ADF\)\s*2018\s*\*.*',
					subMatchers=[
					    SimpleMatcher(startReStr=r'\s\*\s{47}r(?P<program_version>\d+).*')
					],	
					endReStr=r'\s\*{2}.*'),
			    SimpleMatcher(
				    startReStr=r'Geometry.*',
				    sections=['section_system'],
					subMatchers=[
					    SimpleMatcher(         
						    startReStr=r'  Index Symbol       x \(bohr\)       y \(bohr\)       z \(bohr\).*',
						    subMatchers=[
							    SimpleMatcher(
								    startReStr=r'\s*\d+\s+([A-Z][a-z]?)\s+(\d+\.\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)\s*', repeats=True, startReAction=self.save_atoms)
						    ],
							endReStr=r'\s*'
						),
					    SimpleMatcher(
						    startReStr=r'Lattice vectors \(bohr\).*',
						    subMatchers=[
							    SimpleMatcher(
								    startReStr=r'\s*\d+\s+(\d+\.\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)\s*', repeats=True, startReAction=self.save_lattice)
						    ],
							endReStr=r'\s*'
						)	                        
					]	
				),
			    SimpleMatcher(
				    startReStr=r' DENSITY FUNCTIONAL POTENTIAL \(scf\)',
				    sections=['section_method'],
					subMatchers=[
					    #SimpleMatcher(
                        #    sections=['section_XC_functionals'],
						#    startReStr=r'    LDA:                               (?P<XC_functional_name>[a-z,A-Z]*).*'
                        #),
                        SimpleMatcher(
						    startReStr=r'    Gradient Corrections:              ([a-z,A-Z]*)c\s*([a-z,A-Z]*)x.*', startReAction=self.save_functional
                        ),
                        SimpleMatcher(
                            sections=['section_XC_functionals'],
						    startReStr=r'    Meta-GGA:                          (?P<XC_functional_name>[a-z,A-Z]*).*',
                        )                        
					],
                    endReStr=r' DENSITY FUNCTIONAL ENERGY \(post-scf\)'
				),
                SimpleMatcher(
				    sections=['section_single_configuration_calculation'],
				    startReStr=r'Energy \(hartree\)            (?P<energy_total>-\d+\.\d+).*'
				)								
            ]
        )
    def hallo(*args):
        print("hallo")
		
    def save_atoms(self, _, groups):
	    self.atom_positions.append([float(groups[1]), float(groups[2]), float(groups[3])])
	    self.atom_labels.append(groups[0])
        
    def save_lattice(self, _, groups):
	    self.lattice_vectors.append([float(groups[0]), float(groups[1]), float(groups[2])])
        
    def save_functional(self, _, groups):
        if groups != None:
            #print('hallo')
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

    def onClose_section_method(self, backend, *args, **kwargs):
        backend.addValue('electronic_structure_method', 'DFT')
        if self.GGA_functional_name != []:
            print(self.GGA_functional_name)
            backend.openNonOverlappingSection('section_XC_functionals')
            backend.addValue('XC_functional_name', 'GGA_X_' + self.GGA_functional_name[0][1])
            backend.closeNonOverlappingSection('section_XC_functionals')
            backend.openNonOverlappingSection('section_XC_functionals')
            backend.addValue('XC_functional_name', 'GGA_C_' + self.GGA_functional_name[0][0])
            backend.closeNonOverlappingSection('section_XC_functionals')

    def onClose_section_run(self, backend, *args, **kwargs):
        backend.addValue('program_name', 'band')
        backend.addValue('program_basis_set_type', 'slater')

if __name__ == "__main__":
    parser = BANDParser(backend=LocalBackend)
    parser.parse(sys.argv[1])
    #parser.parser_context.super_backend.write_json(sys.stdout)
    print(parser.parser_context.super_backend)
    
#TODO
#total energy Energy (hartree)            
#band gap
#dipole moment  Final bond energy (PBE)
#DOS
#charge
#k-points
#basis set with core treatment
#ADF: Not yet in AMS
