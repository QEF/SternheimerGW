#------------------------------------------------------------------------------
#
# This file is part of the Sternheimer-GW code.
# Parts of this file are taken from the Quantum ESPRESSO software
# P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
#
# Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
# Henry Lambert, Martin Schlipf, and Feliciano Giustino
#
# Sternheimer-GW is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Sternheimer-GW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Sternheimer-GW. If not, see
# http://www.gnu.org/licenses/gpl.html .
#
#------------------------------------------------------------------------------

include ./ENVIRONMENT

default :
	@echo 'To run Sternheimer GW test-suite, type at the shell prompt:'
	@echo ' '
	@echo '  make target'
	@echo ' '
	@echo 'where <target> identifies an action'
	@echo ' run-tests          : run-tests-serial (default)'
	@echo ' run-tests-serial   : run all tests in SERIAL' 
	@echo ' run-tests-parallel : run all tests in PARALLEL (4 MPI)'
	@echo ' compare            : compare last output generated with reference'
	@echo ' clean              : clean output (our/err files) of all tests'
	@echo ' '
	@echo 'For additional advanced commands and settings please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

pseudo :
	@./check_pseudo.sh sgw_

run-tests : run-tests-serial

run-tests-serial : run-tests-sgw-serial

run-tests-parallel : run-tests-sgw-parallel


run-tests-sgw-serial : pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=sgw_all

run-tests-sgw-parallel : pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=sgw_all

run-custom-test : run-custom-test-serial

run-custom-test-serial : pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi

run-custom-test-parallel : pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi


create-reference-sgw :
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=sgw_all make-benchmarks

compare : compare-sgw

compare-sgw :
	${TESTCODE_DIR}/bin/testcode.py --category=sgw_all  --verbose compare 2>&1 | tee out.SGW.`date +%Y%m%d_%H%M%S`

clean:
	@for x in `find sgw_* -name "test*"`; do rm -rf $$x; done
	@for x in `ls -d sgw_*`; do rm -rf $$x/tmp; done

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find sgw_* -name "benchmark*"`; do rm -rf $$x; done
