"""Test the field module"""
import unittest
import tempfile
import sys
import subprocess
import pyopal.objects.parser
import pyopal.objects.field

def run_encapsulated_test(file_name):
    """Run the test in a subfile so that we don't get namespace clashes"""
    pyopal.objects.parser.initialise_from_opal_file(file_name)
    print("Testing")
    value = pyopal.objects.field.get_field_value(0, 1, 0, 0)
    print(value)
    my_return_code = abs(value[3]-2) > 1e-3
    return my_return_code

class FieldTest(unittest.TestCase):
    """Test that we can get fields out"""
    def make_temp(self, a_string):
        """Dump string to a temporary file for use in testing"""
        my_temp = tempfile.NamedTemporaryFile(mode='w+')
        my_temp.write(a_string)
        my_temp.flush()
        return my_temp

    def test_get_field_value(self):
        """Test that we can get out a field value"""
        temp_file = self.make_temp(FieldTest.good_lattice)
        if self.debug:
            temp_stdout = None
            temp_stderr = None
        else:
            temp_stdout = tempfile.TemporaryFile()
            temp_stderr = subprocess.STDOUT
        proc = subprocess.run(["python3", __file__, "run_unit_test_encapsulation", temp_file.name],
                    stdout=temp_stdout, stderr=temp_stderr, check=False)
        if temp_stdout:
            temp_stdout.close() # why does this not use resource allocation?
        self.assertEqual(proc.returncode, 0)

    debug = 0
    command = """
import pyopal.objects.parser
pyopal.objects.parser.initialise_from_opal_file(
"""
    good_lattice = """
///////////////////////////////////////
Title,string="Dummy lattice for testing";
///

OPTION, VERSION=20210100;

field: VERTICALFFAMAGNET,
        B0=2,
        FIELD_INDEX=2,
        MAX_HORIZONTAL_POWER=4,
        END_LENGTH=0.01,
        CENTRE_LENGTH=1000,
        HEIGHT_NEG_EXTENT=4,
        HEIGHT_POS_EXTENT=4,
        WIDTH=2000,
        BB_LENGTH=10;


ringdef: RINGDEFINITION, HARMONIC_NUMBER=1, LAT_RINIT=1, LAT_PHIINIT=0,
         LAT_THETAINIT=0.0, BEAM_PHIINIT=0, BEAM_PRINIT=0,
         BEAM_RINIT=1, SYMMETRY=1, RFFREQ=1, IS_CLOSED=false,
         MIN_R=0.1, MAX_R=2;

l1: Line = (ringdef, field);

Dist1: DISTRIBUTION, TYPE=gauss;

Fs1:FIELDSOLVER, FSTYPE=None, MX=5, MY=5, MT=5,
                 PARFFTX=true, PARFFTY=false, PARFFTT=false,
                 BCFFTX=open, BCFFTY=open, BCFFTT=open,BBOXINCR=2;

beam1: BEAM, PARTICLE=PROTON, pc=1, NPART=1, BCURRENT=1.6e-19, CHARGE=1.0, BFREQ=1;

TRACK, LINE=l1, BEAM=beam1, MAXSTEPS=1, STEPSPERTURN=1;
RUN, METHOD="CYCLOTRON-T", BEAM=beam1, FIELDSOLVER=Fs1, DISTRIBUTION=Dist1;
ENDTRACK;
STOP;
"""

if __name__ == "__main__":
    print(sys.argv, len(sys.argv))
    if len(sys.argv) > 2 and sys.argv[1] == "run_unit_test_encapsulation":
        return_code = run_encapsulated_test(sys.argv[2])
        sys.exit(return_code)
    else:
        FieldTest.debug = 1
        unittest.main()
