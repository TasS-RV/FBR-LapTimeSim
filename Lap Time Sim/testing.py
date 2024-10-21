import numpy as np
from scipy.interpolate import interp1d
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import ezdxf

a = [1,2,3,4,5,6,7,8,9,0]
print(np.take(np.array(a),[3,1,5]))

def colorplot(ax, x, y, c, cmap_name='turbo', vmin=0, vmax=1):
    cmap=mpl.colormaps[cmap_name]
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    for i in range(len(x)-1):
        ax.plot([x[i],x[i+1]],[y[i],y[i+1]],c=cmap(norm(0.5*(c[i]+c[i+1]))))

def read_track(filename):
    import sys
    try:
        doc = ezdxf.readfile(filename)
    except IOError:
        print(f"Not a DXF file or a generic I/O error.")
        sys.exit(1)
    except ezdxf.DXFStructureError:
        print(f"Invalid or corrupted DXF file.")
        sys.exit(2)
    return doc

def find_reverse_segments(entity_list):
    '''finds which segments are actually being traversed backwards, returns a list of reverse=True'''
    start_points=[]
    end_points=[]
    reverse=np.zeros(len(entity_list))
    line_index=0
    for i in range(len(entity_list)):
        if entity_list[i].dxftype()=='ARC':
            center = np.array((entity_list[i].dxfattribs()['center'][0],entity_list[i].dxfattribs()['center'][1]))
            r = entity_list[i].dxfattribs()['radius']
            start_points.append(np.array((center[0]+r*np.cos(entity_list[i].dxfattribs()['start_angle']*np.pi/180),center[1]+r*np.sin(entity_list[i].dxfattribs()['start_angle']*np.pi/180))))
            end_points.append(np.array((center[0]+r*np.cos(entity_list[i].dxfattribs()['end_angle']*np.pi/180),center[1]+r*np.sin(entity_list[i].dxfattribs()['end_angle']*np.pi/180))))
        elif entity_list[i].dxftype()=='LINE':
            line_index=i
            start_points.append(np.array((entity_list[i].dxfattribs()['start'][0],entity_list[i].dxfattribs()['start'][1])))
            end_points.append(np.array((entity_list[i].dxfattribs()['end'][0],entity_list[i].dxfattribs()['end'][1])))
        else:
            print(entity_list[i].dxftype())
            print(entity_list[i].dxfattribs())
    print(line_index)
    #reorder the dxf
    
    i = line_index
    correct_order=[line_index]
    for unused_index in range(len(entity_list)-1):
        delta_start = np.ones(len(entity_list))
        delta_end = np.ones(len((entity_list)))
        if reverse[i]:
            end_point = start_points[i%len(entity_list)]
        else:
            end_point = end_points[i%len(entity_list)]
        for j in range(i+1,i+len(start_points)):
            delta_start[j%len(entity_list)] = (start_points[j%len(entity_list)][0] - end_point[0])**2 +(start_points[j%len(entity_list)][1] - end_point[1])**2
            delta_end[j%len(entity_list)] = (end_points[j%len(entity_list)][0] - end_point[0])**2 +(end_points[j%len(entity_list)][1] - end_point[1])**2
        
        if np.min(delta_start) < 0.01:
            correct_order.append(np.argmin(delta_start))
        elif np.min(delta_end) < 0.01:
            correct_order.append(np.argmin(delta_end))
            reverse[np.argmin(delta_end)] = 1
        else:
            print('some ripe bullshit')
        i = correct_order[-1]
        
    reordered_entity_list = np.take(np.array(entity_list),correct_order)
    reordered_reverse= np.take(reverse,correct_order)

    return reordered_entity_list,  reordered_reverse

fig = plt.figure(figsize=(9,6))
ax = []
ax.append(fig.add_subplot(121))
ax.append(fig.add_subplot(122))

track = read_track('FSA Track.dxf')
entity_list=[e for e in track.modelspace() if e.dxftype()!='MTEXT']
print('{} entities'.format(len(entity_list)))
entity_list, reverse = find_reverse_segments(entity_list)
lengths=[]
for i in range(len(entity_list)):
    entity = entity_list[i]
    if entity.dxftype()=='ARC':
        segment_length = entity.dxf.radius*((entity.dxf.end_angle-entity.dxf.start_angle)%360)*np.pi/180
        lengths.append(segment_length)

        center = np.array((entity.dxfattribs()['center'][0],entity.dxfattribs()['center'][1]))
        r = entity.dxfattribs()['radius']
        if entity.dxfattribs()['start_angle']<entity.dxfattribs()['end_angle']:
            theta = np.linspace(entity.dxfattribs()['start_angle'],entity.dxfattribs()['end_angle'])*np.pi/180
        else:
            theta = np.linspace(entity.dxfattribs()['start_angle'],entity.dxfattribs()['end_angle']+360)*np.pi/180
        if reverse[i]:
            theta = np.flip(theta)

        xmap = center[0]+r*np.cos(theta)
        ymap = center[1]+r*np.sin(theta)

    elif entity.dxftype()=='LINE':
        segment_length=np.linalg.norm([entity.dxf.end[0]-entity.dxf.start[0],entity.dxf.end[1]-entity.dxf.start[1]])
        lengths.append(segment_length)

        start = np.array((entity.dxfattribs()['start'][0],entity.dxfattribs()['start'][1]))
        end = np.array((entity.dxfattribs()['end'][0],entity.dxfattribs()['end'][1]))
        xmap = np.linspace(start[0],end[0])
        ymap = np.linspace(start[1],end[1])
    else:
        print(entity.dxftype())
    colorplot(ax[1],xmap,ymap,i*np.ones(len(xmap)), vmin=0, vmax=43)
    start_x=np.sum(lengths)-segment_length
    end_x = np.sum(lengths)
    ax[0].plot([start_x,end_x],[i,i+1],color=mpl.colormaps['turbo'](mpl.colors.Normalize(vmin=0,vmax=43)(i)))

ax[0].grid()
ax[1].set(aspect='equal')
#cax = plt.axes([0.85, 0.1, 0.075, 0.8])
sm = plt.cm.ScalarMappable(cmap=mpl.colormaps['jet'], norm = mpl.colors.Normalize(vmin=0,vmax=43))
sm.set_array([])
plt.colorbar(sm, ticks=np.linspace(0, 43, 10), ax = ax[0])
plt.show()