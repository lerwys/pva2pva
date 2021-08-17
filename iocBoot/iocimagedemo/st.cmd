#!../../bin/linux-x86_64/softIocPVA

dbLoadRecords("image.db","N=TST:image1")
dbLoadRecords("ndain.db","N=TST:image2")
dbLoadRecords("ndaroi.db","N=TST:image3")
dbLoadRecords("table.db","N=TST:table1")

iocInit()
