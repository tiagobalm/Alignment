import requests
import re
from datetime import datetime
from alignment import needleman_wunsch

# EBI URLS
EBI_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/"
EBI_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"
REGEX = "(?<=\d )[AGVLYETDFQHICSMKPRNW\-]{1,50}"
FILE1 = "HumanHemoglobin"
FILE2 = "AtelesHemoglobin"

# Request data
data = {"email": "up201305665@fe.up.pt",
        "matrix": "EBLOSUM62",
        "gapopen": "10.0",
        "gapext": "0.5",
        "endweight": "false",
        "endopen": "10.0",
        "endextend": "0.5",
        "format": "pair",
        "stype": "protein",
        "asequence": open("./samples/"+FILE1, "r").read(),
        "bsequence": open("./samples/"+FILE2, "r").read()
        }


def get_alignment_from_ebi():

    print("[{0}] SENDING REQUEST TO {1}".format(datetime.now(),EBI_URL))
    r = requests.post(EBI_URL, data)
    print ("JOB ID = {}".format(r.text))

    result = requests.get(EBI_RESULT_URL + "{jobID}/aln".format(jobID=r.text))

    while result.status_code is not 200:
        print("[{0}] NOT READY".format(datetime.now()))
        result = requests.get(EBI_RESULT_URL + "{jobID}/aln".format(jobID=r.text))

    result = re.sub(re.compile("#.*?\n"), "", result.text)

    print("[{0}] WEBSERVICE RESULT {1}".format(datetime.now(), result))

    result = result.split('\n')
    alignmentEBI = result[3::4]
    alignmentEBI = ''.join(alignmentEBI).replace(" ", "")

    asequence = open("./samples/" + FILE1, "r").read()
    bsequence = open("./samples/" + FILE2, "r").read()

    asequence = asequence.split('\n')[1::]
    bsequence = bsequence.split('\n')[1::]

    asequence = ''.join(asequence)
    bsequence = ''.join(bsequence)

    print("EBI RESULT: ")
    print(asequence)
    print(alignmentEBI)
    print(bsequence)
    print()

    print('-' * 100)
    print()
    print("LOCAL RESULT: ")
    needleman_wunsch(asequence, bsequence, -10, False)


get_alignment_from_ebi()
