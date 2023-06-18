import kp

k = kp.kp()

#Ask for data file name
datafile = raw_input(prompt="Name of data file (usually batch_t450.dat or single_t450.dat): ")

#load data
k.loadDataFromFile(datafile)

#set xlabel
k.xlabel = r'density/cm$^{-3}$'

#plot type
k.plog = "loglog"

#plot
k.ylabel = "fraction"
k.ymax = 3e0
k.multiplotLayout("21",sharex=True)
k.multiplotAdd(columns=["ntot","QH","QH2"])

k.ymin = 1e2
k.ylabel = "T/K"
k.multiplotAdd(columns=["ntot","Tgas"])
k.multiplotShow(outputFileName="plot.png")

