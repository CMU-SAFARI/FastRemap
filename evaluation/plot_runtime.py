import glob 
import matplotlib.pyplot as plt 
import numpy as np 

human = ['hg16_hg17', 'hg16_hg18', 'hg16_hg19', 'hg16_hg38', 'hg17_hg18', 'hg17_hg19', 'hg17_hg38', 'hg18_hg19', 'hg18_hg38', 'hg19_hg38'] 
ce = ['ce2_ce4', 'ce4_ce6', 'ce4_ce10', 'ce6_ce10', 'ce6_ce11', 'ce10_ce11'] 
yeast = ['sacCer1_sacCer3', 'sacCer2_sacCer3'] 

def compare_arrays(x, y, y2): 
    for idx, el in enumerate(y): 
        print(x[idx]) 
        print(str(y[idx]) + " : " + str(y2[idx])  + "  :  " + str(float(y[idx]) / float(y2[idx]))) 
        

def wall_clock_to_seconds(wallclock):
    time = 0 
    for idx, el in enumerate(reversed(wallclock.split(":"))): 
        if (idx == 0):
            time += float(el) 
        if (idx == 1): 
            time += float(el) * 60 
        if (idx == 2): 
            time += float(el) * 60 * 60 
        if (idx == 3): 
            print("wall_clock_to_seconds broken") 
            exit() 

    return time 

crossmap  = dict() 
fastremap = dict() 

for fn in glob.glob("./crossmap/*.time"): 
    pair = fn.split("/")[-1].split(".")[0] 
    if pair not in crossmap: 
        crossmap[pair] = [0, 0, 0, 0] 
    for line in open(fn): 
        if "User time" in line: 
            crossmap[pair][0] = float(line.split()[-1]) 
        elif "System time" in line: 
            crossmap[pair][1] = float(line.split()[-1])
        elif "wall clock)" in line: 
            crossmap[pair][2] = wall_clock_to_seconds(line.split()[-1]) 
        elif "Maximum resident set size" in line: 
            crossmap[pair][3] = float(line.split()[-1]) 

for fn in glob.glob("./FastRemap/*.time"): 
    pair = fn.split("/")[-1].split(".")[0] 
    if pair not in fastremap: 
        fastremap[pair] = [0, 0, 0, 0] 
    for line in open(fn): 
        if "User time" in line: 
            fastremap[pair][0] = float(line.split()[-1]) 
        elif "System time" in line: 
            fastremap[pair][1] = float(line.split()[-1])
        elif "wall clock)" in line: 
            fastremap[pair][2] = wall_clock_to_seconds(line.split()[-1]) 
        elif "Maximum resident set size" in line: 
            fastremap[pair][3] = float(line.split()[-1]) 
    
    print(fn) 

print(crossmap) 
print(fastremap) 

width = 0.35

x = [] 
y_sys_user = [] 
y_wallclock = [] 
y_mem = [] 
y2_sys_user = [] 
y2_wallclock = [] 
y2_mem = [] 


fig_h1,ax_h1 = plt.subplots(figsize=[10,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
fig_h2,ax_h2 = plt.subplots(figsize=[10,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
fig_h3,ax_h3 = plt.subplots(figsize=[10,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
# human 
for el in human: 
    if el not in crossmap or el not in fastremap: 
        continue 
    x.append(el) 
    y_sys_user.append(crossmap[el][0] + crossmap[el][1]) 
    y_wallclock.append(crossmap[el][2]) 
    y_mem.append(crossmap[el][3]) 

    y2_sys_user.append(fastremap[el][0] + fastremap[el][1]) 
    y2_wallclock.append(fastremap[el][2]) 
    y2_mem.append(fastremap[el][3]) 

x_axis = np.arange(len(x)) 
ax_h1.bar(x_axis - width/2, y_sys_user, width, label="Crossmap") 
ax_h1.bar(x_axis + width/2, y2_sys_user, width, label="FastRemap") 
ax_h1.set_xticks(x_axis) 
ax_h1.set_xticklabels(x) 
ax_h1.set_title("h sys_user time") 
#ax_h1.legend() 
fig_h1.savefig("runtime/h_sys_user_time.pdf") 
print("H SYS USER TIME") 
compare_arrays(x, y_sys_user, y2_sys_user)

ax_h2.bar(x_axis - width/2, y_wallclock, width, label="Crossmap") 
ax_h2.bar(x_axis + width/2, y2_wallclock, width, label="FastRemap") 
ax_h2.set_xticks(x_axis) 
ax_h2.set_xticklabels(x) 
ax_h2.set_title("h wallclock time") 
#ax_h2.legend() 
fig_h2.savefig("runtime/h_wallclock_time.pdf") 
print("H WALLCLOCK TIME") 
compare_arrays(x, y_wallclock, y2_wallclock)

ax_h3.bar(x_axis - width/2, y_mem, width, label="Crossmap") 
ax_h3.bar(x_axis + width/2, y2_mem, width, label="FastRemap") 
ax_h3.set_xticks(x_axis) 
ax_h3.set_xticklabels(x) 
ax_h3.set_title("h mem") 
#ax_h3.legend() 
fig_h3.savefig("runtime/h_mem.pdf") 
print("H MEM") 
compare_arrays(x, y2_mem, y_mem)

x.clear(); y_sys_user.clear(); y_wallclock.clear(); y_mem.clear(); y2_sys_user.clear(); y2_wallclock.clear(); y2_mem.clear() 


fig_c1,ax_c1 = plt.subplots(figsize=[6,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
fig_c2,ax_c2 = plt.subplots(figsize=[6,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
fig_c3,ax_c3 = plt.subplots(figsize=[6,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
# ce 
for el in ce: 
    if el not in crossmap or el not in fastremap: 
        continue 
    x.append(el) 
    y_sys_user.append(crossmap[el][0] + crossmap[el][1]) 
    y_wallclock.append(crossmap[el][2]) 
    y_mem.append(crossmap[el][3]) 

    y2_sys_user.append(fastremap[el][0] + fastremap[el][1]) 
    y2_wallclock.append(fastremap[el][2]) 
    y2_mem.append(fastremap[el][3]) 

x_axis = np.arange(len(x)) 
ax_c1.bar(x_axis - width/2, y_sys_user, width, label="Crossmap") 
ax_c1.bar(x_axis + width/2, y2_sys_user, width, label="FastRemap") 
ax_c1.set_xticks(x_axis) 
ax_c1.set_xticklabels(x) 
#ax_c1.legend() 
ax_c1.set_title("ce sys_user time") 
fig_c1.savefig("runtime/c_sys_user_time.pdf") 
print("CE SYS USER TIME") 
compare_arrays(x, y_sys_user, y2_sys_user)


ax_c2.bar(x_axis - width/2, y_wallclock, width, label="Crossmap") 
ax_c2.bar(x_axis + width/2, y2_wallclock, width, label="FastRemap") 
ax_c2.set_xticks(x_axis) 
ax_c2.set_xticklabels(x) 
#ax_c2.legend() 
ax_c2.set_title("ce wallclock time") 
fig_c2.savefig("runtime/c_wallclock_time.pdf") 
print("CE WALLCLOCK TIME") 
compare_arrays(x, y_wallclock, y2_wallclock)

ax_c3.bar(x_axis - width/2, y_mem, width, label="Crossmap") 
ax_c3.bar(x_axis + width/2, y2_mem, width, label="FastRemap") 
ax_c3.set_xticks(x_axis) 
ax_c3.set_xticklabels(x) 
#ax_c3.legend() 
ax_c3.set_title("ce mem") 
fig_c3.savefig("runtime/c_mem.pdf") 
print("CE MEM") 
compare_arrays(x, y2_mem, y_mem)



x.clear(); y_sys_user.clear(); y_wallclock.clear(); y_mem.clear(); y2_sys_user.clear(); y2_wallclock.clear(); y2_mem.clear() 

fig_s1,ax_s1 = plt.subplots(figsize=[2,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
fig_s2,ax_s2 = plt.subplots(figsize=[2,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
fig_s3,ax_s3 = plt.subplots(figsize=[2,5]) 
plt.grid(linestyle='--', linewidth=0.5, axis='y') 
# saccer
for el in yeast: 
    if el not in crossmap or el not in fastremap: 
        continue 
    x.append(el) 
    y_sys_user.append(crossmap[el][0] + crossmap[el][1]) 
    y_wallclock.append(crossmap[el][2]) 
    y_mem.append(crossmap[el][3]) 

    y2_sys_user.append(fastremap[el][0] + fastremap[el][1]) 
    y2_wallclock.append(fastremap[el][2]) 
    y2_mem.append(fastremap[el][3]) 

x_axis = np.arange(len(x)) 
ax_s1.bar(x_axis - width/2, y_sys_user, width, label="Crossmap") 
ax_s1.bar(x_axis + width/2, y2_sys_user, width, label="FastRemap") 
ax_s1.set_xticks(x_axis) 
ax_s1.set_xticklabels(x) 
#ax_s1.legend() 
ax_s1.set_title("yeast sys_user time") 
fig_s1.savefig("runtime/s_sys_user_time.pdf") 
print("SACCER SYS USER TIME") 
compare_arrays(x, y_sys_user, y2_sys_user)

ax_s2.bar(x_axis - width/2, y_wallclock, width, label="Crossmap") 
ax_s2.bar(x_axis + width/2, y2_wallclock, width, label="FastRemap") 
ax_s2.set_xticks(x_axis) 
ax_s2.set_xticklabels(x) 
#ax_s2.legend() 
ax_s2.set_title("yeast wallclock time") 
fig_s2.savefig("runtime/s_wallclock_time.pdf") 
print("SACCER WALLCLOCK TIME") 
compare_arrays(x, y_wallclock, y2_wallclock)

ax_s3.bar(x_axis - width/2, y_mem, width, label="Crossmap") 
ax_s3.bar(x_axis + width/2, y2_mem, width, label="FastRemap") 
ax_s3.set_xticks(x_axis) 
ax_s3.set_xticklabels(x) 
#ax_s3.legend() 
ax_s3.set_title("yeast mem") 
fig_s3.savefig("runtime/s_mem.pdf") 
print("SACCER MEM") 
compare_arrays(x, y2_mem, y_mem)

plt.show() 

