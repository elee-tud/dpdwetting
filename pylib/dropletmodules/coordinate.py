#!/home/elee/programs/Python-3.7.10/bin/python

def write_gro(filename, title, topol, box, coord, vel):
    fpout=open(filename, "w");
    nbeads=len(topol)
    fpout.write("%s\n"%title)
    fpout.write("%d\n"%nbeads)

    for i in range(nbeads):
        fpout.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"%(topol[i][0]%100000, topol[i][1], topol[i][2],  topol[i][3]%100000, coord[i][0], coord[i][1], coord[i][2], vel[i][0], vel[i][1], vel[i][2]))
    fpout.write("%10.5f%10.5f%10.5f\n"%(box[0], box[1], box[2]))
    fpout.close()


def read_gro(filename):
    topol=[]
    box=[]
    coord=[]
    vel=[]
    fpin=open(filename, "r")
    data=fpin.readlines()
    title=str(data[0][:-1]).rstrip("\n")
    for line in data[2:-1]:
        molnum=int(line[0:5])
        molname=str(line[5:10])
        attype=str(line[10:15])
        atnum=int(line[15:20])
        x=float(line[20:28])
        y=float(line[28:36])
        z=float(line[36:44])
        vx=float(line[44:52])
        vy=float(line[52:60])
        vz=float(line[60:68])
        topol.append([molnum, molname, attype, atnum])
        coord.append([x,y,z])
        vel.append([vx,vy,vz])
    box.append(float(data[-1][0:10]))
    box.append(float(data[-1][10:20]))
    box.append(float(data[-1][20:30]))



    return title, topol, box, coord, vel

def read_gro_att(fpin, skip=False):
    topol=[]
    box=[]
    coord=[]
    vel=[]
    title=str(fpin.readline()).rstrip("\n")
    natoms=int(fpin.readline())
    for i in range(natoms):
        line=fpin.readline()
        if not skip:
            molnum=int(line[0:5])
            molname=str(line[5:10])
            attype=str(line[10:15])
            atnum=int(line[15:20])
            x=float(line[20:28])
            y=float(line[28:36])
            z=float(line[36:44])
            vx=float(line[44:52])
            vy=float(line[52:60])
            vz=float(line[60:68])
            topol.append([molnum, molname, attype, atnum])
            coord.append([x,y,z])
            vel.append([vx,vy,vz])
    line=fpin.readline()
    if not skip:
        box.append(float(line[0:10]))
        box.append(float(line[10:20]))
        box.append(float(line[20:30]))

    return title, topol, box, coord, vel

