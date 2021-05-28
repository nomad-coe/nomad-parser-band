#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from bandparser import BandParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return BandParser()


def test_scf(parser):
    archive = EntryArchive()
    parser.parse('tests/data/phenylrSmall-metagga.out', archive, None)

    sec_run = archive.section_run[0]
    assert sec_run.program_version == '70630 2018-11-24'
    assert sec_run.time_run_date_start.magnitude == 1550053118.0

    sec_method = archive.section_run[0].section_method[0]
    assert sec_method.number_of_spin_channels == 2
    assert sec_method.section_XC_functionals[0].XC_functional_name == 'MGGA_XC_TPSS'

    sec_system = archive.section_run[0].section_system[0]
    assert sec_system.atom_positions[2][1].magnitude == approx(2.31103866e-10)
    assert sec_system.lattice_vectors[0][0].magnitude == approx(8.72113987e-10)
    assert sec_system.configuration_periodic_dimensions == [True, True, False]

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert sec_scc.energy_total.magnitude == approx(-1.11415974e-17)
    assert sec_scc.electronic_kinetic_energy.magnitude == approx(1.11552116e-17)
    assert sec_scc.energy_XC.magnitude == approx(-1.01588269e-17)
    assert sec_scc.energy_electrostatic.magnitude == approx(-8.69562185e-18)

    sec_scfs = sec_scc.section_scf_iteration
    assert len(sec_scfs) == 20
    assert sec_scfs[11].energy_change_scf_iteration.magnitude == approx(2.1449944e-21)


def test_geometry_optimization(parser):
    archive = EntryArchive()
    parser.parse('tests/data/phenylrSmall-geoopt.out', archive, None)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 22
    assert sec_sccs[17].energy_total.magnitude == approx(-1.10701399e-17)
    assert len(sec_sccs[6].section_scf_iteration) == 11
    assert sec_sccs[10].atom_forces[5][1].magnitude == approx(3.11423748e-11)

    sec_systems = archive.section_run[0].section_system
    assert sec_systems[3].atom_positions[10][1].magnitude == approx(5.81904209e-10)
    assert sec_systems[9].lattice_vectors[1][1].magnitude == approx(7.56032016e-10)


def test_dos(parser):
    archive = EntryArchive()
    parser.parse('tests/data/NiO-dos.out', archive, None)

    sec_dos = archive.section_run[0].section_single_configuration_calculation[0].dos_electronic[0]
    assert np.shape(sec_dos.dos_total[1].dos_values) == (158,)
    assert sec_dos.dos_energies[78].magnitude == approx(-8.6717613e-20)
    assert sec_dos.dos_total[1].dos_values[19] == approx(4.66493971e+14)

    archive = EntryArchive()
    parser.parse('tests/data/NiO-dos-restricted.out', archive, None)
    sec_dos = archive.section_run[0].section_single_configuration_calculation[0].dos_electronic[0]
    assert np.shape(sec_dos.dos_total[0].dos_values) == (154,)
