from tkinter import *
from tkinter import ttk
import tkinter.scrolledtext as tkst
from alignmentEBI import get_alignment_from_ebi
from alignment import align_locally
from gapPenaltyAlignment import gap_penalty_align


def get_result(event):
    status.set("PROCESSING...")
    Tk.update(root)

    data = {"email": "up201305665@fe.up.pt",
            "matrix": "E"+defaultMatrix.get(),
            "gapopen": gapOpenVariable.get(),
            "gapext": gapNextVariable.get(),
            "endweight": "true" if endWeightVar.get() else "false",
            "endopen": endOpenVariable.get(),
            "endextend": endExtendVariable.get(),
            "format": defaultFormat.get().lower(),
            "stype": defaultStyle.get().lower(),
            "asequence": asequenceText.get(1.0, END),
            "bsequence": bsequenceText.get(1.0, END)
            }

    ebi_result = '\n'.join(get_alignment_from_ebi(data))
    ebiResult.delete(1.0, END)
    ebiResult.insert(INSERT, ebi_result)

    local_result = ""

    if radioButtonVar.get() == "linear":
        local_results = align_locally(asequenceText.get(1.0, END), bsequenceText.get(1.0, END),
                                     float(gapOpenVariable.get()) * -1,
                                      True if defaultStyle.get().lower() == "dna" else False)
    else:
        local_results = gap_penalty_align(asequenceText.get(1.0, END), bsequenceText.get(1.0, END),
                                      float(gapOpenVariable.get()) * -1,
                                          float(gapNextVariable.get()) * -1,
                                          True if defaultStyle.get().lower() == "dna" else False)

    for result in local_results:
        local_result += '\n'.join(result)
        local_result += '\n\n'

    localResult.delete(1.0, END)
    localResult.insert(INSERT, local_result)

    status.set("FINISHED")

    return


def disable_input():
    gapNextText.configure(state='disabled')
    endExtendText.configure(state='disabled')
    endWeightButton.configure(state='disabled')
    endOpenText.configure(state='disabled')
    return


def enable_input():
    gapNextText.configure(state='normal')
    endWeightButton.configure(state='normal')

    if endWeightVar.get():
        endExtendText.configure(state='normal')
        endOpenText.configure(state='normal')

    return


def toggle_end_options():
    if endWeightVar.get():
        endExtendText.configure(state='normal')
        endOpenText.configure(state='normal')
    else:
        endExtendText.configure(state='disabled')
        endOpenText.configure(state='disabled')
    return


root = Tk()

root.geometry("1280x720+300+150")
root.resizable(width=False, height=False)
root.title("Pairwise Alignment")

mainFrame = Frame(root)

sequenceFrame = Frame(mainFrame)

asequenceLabel = Label(sequenceFrame, text="First Sequence (Fasta)")
asequenceLabel.grid(row=0, column=0, sticky=W)
asequenceText = tkst.ScrolledText(sequenceFrame, width=90, height=8, wrap='word')
asequenceText.grid(row=1, column=0, sticky=W)

bsequenceLabel = Label(sequenceFrame, text="Second Sequence (Fasta)")
bsequenceLabel.grid(row=2, column=0, sticky=W)
bsequenceText = tkst.ScrolledText(sequenceFrame, width=90, height=8, wrap='word')
bsequenceText.grid(row=3, column=0, sticky=W)

optionsFrame = Frame(mainFrame)

radioButtonVar = StringVar()
radioButtonVar.set("linear")

linearCostLabel = Label(optionsFrame, text="Linear Cost")
linearCostLabel.grid(row=0, column=0, sticky=W, pady=(40, 10))
linearCostButton = Radiobutton(optionsFrame, variable=radioButtonVar, value="linear", command=lambda: disable_input())
linearCostButton.grid(row=0, column=1, sticky=W, pady=(40, 10))

gapPenaltyLabel = Label(optionsFrame, text="Gap Penalty")
gapPenaltyLabel.grid(row=0, column=2, sticky=W, pady=(40, 10))
gapPenaltyButton = Radiobutton(optionsFrame, variable=radioButtonVar, value="gap", command=lambda: enable_input())
gapPenaltyButton.grid(row=0, column=3, sticky=W, pady=(40, 10))

matrixLabel = Label(optionsFrame, text="Matrix")
matrixLabel.grid(row=1, column=0, sticky=W, pady=(0, 10))

defaultMatrix = StringVar(optionsFrame)
defaultMatrix.set("BLOSUM62")
matrixDropDown = OptionMenu(optionsFrame, defaultMatrix, "BLOSUM20")
matrixDropDown.grid(row=1, column=1, sticky=(W, E), pady=(0, 10))

gapOpenLabel = Label(optionsFrame, text="GAPOPEN")
gapOpenLabel.grid(row=2, column=0, sticky=W, pady=(0, 10))
gapOpenVariable = StringVar()
gapOpenVariable.set("10.0")
gapOpenText = Entry(optionsFrame, textvariable=gapOpenVariable)
gapOpenText.grid(row=2, column=1, sticky=W, pady=(0, 10))

gapNextLabel = Label(optionsFrame, text="GAPNEXT")
gapNextLabel.grid(row=2, column=2, sticky=W, pady=(0, 10), padx=(10, 0))
gapNextVariable = StringVar()
gapNextVariable.set("0.5")
gapNextText = Entry(optionsFrame, textvariable=gapNextVariable)
gapNextText.grid(row=2, column=3, sticky=W, pady=(0, 10))

endWeight = Label(optionsFrame, text="ENDWEIGHT")
endWeight.grid(row=3, column=0, sticky=W, pady=(0, 10))
endWeightVar = IntVar()
endWeightButton = Checkbutton(optionsFrame, variable=endWeightVar, command=toggle_end_options)
endWeightButton.grid(row=3, column=1, sticky=W, pady=(0, 10))

endOpenLabel = Label(optionsFrame, text="ENDOPEN")
endOpenLabel.grid(row=4, column=0, sticky=W, pady=(0, 10))
endOpenVariable = StringVar()
endOpenVariable.set("10.0")
endOpenText = Entry(optionsFrame, textvariable=endOpenVariable)
endOpenText.grid(row=4, column=1, sticky=W, pady=(0, 10))
endOpenText.configure(state='disabled')

endExtendLabel = Label(optionsFrame, text="ENDEXTEND")
endExtendLabel.grid(row=4, column=2, sticky=W, pady=(0, 10), padx=(10, 0))
endExtendVariable = StringVar()
endExtendVariable.set("0.5")
endExtendText = Entry(optionsFrame, textvariable=endExtendVariable)
endExtendText.grid(row=4, column=3, sticky=W, pady=(0, 10))
endExtendText.configure(state='disabled')

formatLabel = Label(optionsFrame, text="Format")
formatLabel.grid(row=5, column=0, sticky=W, pady=(0, 10))

defaultFormat = StringVar(optionsFrame)
defaultFormat.set("Pair")
formatDropDown = OptionMenu(optionsFrame, defaultFormat, "Pair")
formatDropDown.grid(row=5, column=1, sticky=(W, E), pady=(0, 10))

styleLabel = Label(optionsFrame, text="Style")
styleLabel.grid(row=5, column=2, sticky=W, pady=(0, 10), padx=(10, 0))

defaultStyle = StringVar(optionsFrame)
defaultStyle.set("PROTEIN")
styleDropDown = OptionMenu(optionsFrame, defaultStyle, "PROTEIN", "DNA")
styleDropDown.grid(row=5, column=3, sticky=(W, E), pady=(0, 10))

submitButton = Button(optionsFrame, text="Submit")
submitButton.grid(row=6, column=0, sticky=(W, E), columnspan=4, pady=(30, 0))
submitButton.bind("<Button-1>", get_result)

resultsFrame = Frame(mainFrame)

ebiResultLabel = Label(resultsFrame, text="EBI RESULT")
ebiResultLabel.grid(row=0, column=0, sticky=W, pady=(0, 5))
status = StringVar()
status.set("NOT SUBMITTED")
statusLabel = Label(resultsFrame, textvariable=status)
statusLabel.grid(row=0, column=1, sticky=E, pady=(0, 5))

ebiResult = tkst.ScrolledText(resultsFrame, width=150, height=8, wrap='word')
ebiResult.grid(row=1, column=0, sticky=W, pady=(0, 20))

localResultLabel = Label(resultsFrame, text="LOCAL RESULT")
localResultLabel.grid(row=2, column=0, sticky=W, pady=(0, 5))
localResult = tkst.ScrolledText(resultsFrame, width=150, height=8, wrap='word')
localResult.grid(row=3, column=0, sticky=W)

sequenceFrame.grid(row=0, column=0, sticky=W, padx=(0, 50))
optionsFrame.grid(row=0, column=1, sticky=N, padx=(0, 50))
resultsFrame.grid(row=1, column=0, sticky=S, columnspan=2, pady=(50, 0))
mainFrame.grid(row=0, column=0, sticky=(N, S, W, E), padx=20, pady=20)

root.mainloop()