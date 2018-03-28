import requests
import re
from datetime import datetime
from alignment import needleman_wunsch

# EBI URLS
EBI_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/"
EBI_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"
REGEX = "(?<=\d )[AGVLYETDFQHICSMKPRNW\-]{1,50}"
FILE1 = "HBA_CTEGU"
FILE2 = "HBA_ATEGE"

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


def parse_ebi_result(result):

    result = re.sub(re.compile("#.*?\n"), "", result)
    p = re.compile(REGEX)
    m = p.findall(result)
    asequence = m[0::2]
    bsequence = m[1::2]
    asequence = ''.join(asequence)
    bsequence = ''.join(bsequence)

    alignment = result.split('\n');
    alignment = alignment[3::4]

    for x in range(0, len(alignment)):
        alignment[x] = alignment[x][21::]

    alignment = ''.join(alignment)

    print("EBI RESULT: ")
    print(asequence)
    print(alignment)
    print(bsequence)


def parse_files_for_local_alignment():
    asequence = open("./samples/" + FILE1, "r").read()
    bsequence = open("./samples/" + FILE2, "r").read()

    asequence = asequence.split('\n')[1::]
    bsequence = bsequence.split('\n')[1::]

    asequence = ''.join(asequence)
    bsequence = ''.join(bsequence)

    return asequence, bsequence


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

    parse_ebi_result(result)

    print()

    asequence, bsequence = parse_files_for_local_alignment()

    print('-' * 100)
    print()
    print("LOCAL RESULT: ")
    needleman_wunsch(asequence, bsequence, -10, False)


get_alignment_from_ebi()
