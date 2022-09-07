import os
import sys
def makeproject(PROJECT,SAMPLE):
    os.system(f"mkdir projects/{PROJECT}/main_result/ -p")

    
    os.system(f"mkdir projects/{PROJECT}/sample_{SAMPLE}/snippy -p")
    os.system(f"mkdir projects/{PROJECT}/sample_{SAMPLE}/main_result -p")

    os.system(f"cp samples/{SAMPLE}/snippy/* -T projects/{PROJECT}/sample_{SAMPLE}/snippy/")

makeproject(sys.argv[1], sys.argv[2])