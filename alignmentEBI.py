import requests
import re
from datetime import datetime

# EBI URLS
EBI_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/run/"
EBI_RESULT_URL = "http://www.ebi.ac.uk/Tools/services/rest/emboss_needle/result/"
REGEX = "(?<=\d )[AGVLYETDFQHICSMKPRNW\-]{1,50}"
FILE1 = "HBA_CTEGU"
FILE2 = "HBA_ATEGE"

# Request data example
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
    ebi_result = []

    result = re.sub(re.compile("#.*?\n"), "", result)
    p = re.compile(REGEX)
    m = p.findall(result)
    asequence = m[0::2]
    bsequence = m[1::2]
    asequence = ''.join(asequence)
    bsequence = ''.join(bsequence)

    alignment = result.split('\n')
    alignment = alignment[3::4]

    for x in range(0, len(alignment)):
        alignment[x] = alignment[x][21::]

    alignment = ''.join(alignment)

    ebi_result.append(asequence)
    ebi_result.append(alignment)
    ebi_result.append(bsequence)

    return ebi_result


def get_alignment_from_ebi(request_data):

    print("[{0}] SENDING REQUEST TO {1}".format(datetime.now(), EBI_URL))
    r = requests.post(EBI_URL, request_data)
    print("JOB ID = {}".format(r.text))

    result = requests.get(EBI_RESULT_URL + "{jobID}/aln".format(jobID=r.text))

    while result.status_code is not 200:
        print("[{0}] NOT READY".format(datetime.now()))
        result = requests.get(EBI_RESULT_URL + "{jobID}/aln".format(jobID=r.text))

    result = re.sub(re.compile("#.*?\n"), "", result.text)

    print("[{0}] WEBSERVICE RESULT {1}".format(datetime.now(), result))

    return parse_ebi_result(result)
