#!/usr/bin/env python2
# *- coding: utf-8 -*-
# vim: set expandtab tabstop=4 shiftwidth=4 softtabstop=4 textwidth=79:

from subprocess import Popen, PIPE
import os


class CircosController:
    circos_chr_name = 'hs'# make this argumentable
    smoothing_offset = 15000
    smoothing_prec = 3
    
    def __init__(self, sid, data, main_config_file, coordinate_config_file, data_file):
        self.logger = logging.getLogger(self.__class__.__name__)
        
        self.sid = sid
        self.data = data
        #self.main_config_file = main_config_file
        self.coordinate_config_file = coordinate_config_file
        self.data_file = data_file
    
    def draw_network(self, png_file, svg_file):
        self.write_configs()
    
    def write_configs(self):
        # Chooses the chromosomes
        coordinates = self.estimate_coordinates()
        self.write_coordinate_config(coordinates)#"tmp/select-coordinates.conf")
        #self.write_basic_config()
        
        self.write_data()
        
        self.run()
    
    def write_coordinate_config(self, coordinates):
        """Generates a config file similar to this:
chromosomes_display_default = no
chromosomes                 = hs21[a]:0-39.7;hs21[b]:39.7-39.9;hs21[c]:39.9-42.7;hs21[d]:42.7-42.9;hs21[e]:42.9-50
#chromosomes_reverse         = /hs[]/
#chromosomes_scale           = hs1=0.5r,/hs[234]/=1.0rn
chromosomes_scale           = a:1.0;b:25.0;c:1.0;d:25.0;e:1.0
chromosomes_radius          = a=0.4r,b=0.99r,c=0.4r,d=0.99r,e=0.4r
        """
        
        scales = {'small': '1.0', 'large': '35.0'}
        radius = {'small': '0.55r', 'large': '0.99r'}
        
        fh = open(self.coordinate_config_file, "w")
        fh.write("chromosomes_display_default = no\n")
        fh.write("chromosomes = "+";".join([ c[0]+"["+c[4]+"]:"+(("%."+str(self.smoothing_prec)+"f") % c[1])+"-"+(("%."+str(self.smoothing_prec)+"f") %c[2]) for c in coordinates])+"\n")
        fh.write("chromosome_scale = "+";".join( [str(c[4])+":"+scales[c[3]] for c in coordinates])+"\n")
        fh.write("chromosome_radius = "+",".join( [str(c[4])+"="+radius[c[3]] for c in coordinates])+"\n")
        
        fh.close()

    #def write_basic_config(self):
    #    fh = open(self.main_config_file, "w")
    #    fh.write("")
    #    fh.close()
    
    def estimate_coordinates(self):
        """Returns chunks with a certain offset (usually 10kb)
        """
        
        idx = {}
        for dp in self.data:
            for i in [0,1]:
                _chr = dp[i]._target.position._chr.replace('chr',self.circos_chr_name)
                if not idx.has_key(_chr):
                    idx[_chr] = {}
                
                idx[_chr][dp[i]._target.position.pos] = True
        
        vec = []
        for _chr in idx.keys():
            previous = None
            vec_large = []
            
            for pos in sorted(idx[_chr].keys()):
                chunk = [pos-self.smoothing_offset,pos+self.smoothing_offset]
                chunk[0] = float(chunk[0]) / float(pow(10,6))# a million -> Mb
                chunk[1] = float(chunk[1]) / float(pow(10,6))# a million -> Mb
                chunk = [round(chunk[0],self.smoothing_prec), round(chunk[1],self.smoothing_prec)]
                
                # Look what to do with previous chunk, if there is any
                if previous != None:
                    # if overlap, extend previous chunk
                    if chunk[0] <= previous[1]:
                        chunk[0] = previous[0]
                    else:
                        vec_large.append( (_chr, max(0.0, previous[0]), previous[1], 'large') )
                
                previous = chunk
            
            vec_large.append( (_chr, max(0.0,previous[0]), previous[1], 'large') )
            
            i = 0
            if vec_large[0][0] > 0.0:
                vec.append( (_chr, 0.0, vec_large[0][1], 'small', _chr+"_"+str(i)) )
                i += 1
            
            for k in range(len(vec_large)-1):
                vec.append( (_chr, vec_large[k][1], vec_large[k][2], vec_large[k][3], _chr+"_"+str(i)))
                i += 1
                vec.append( (_chr, vec_large[k][2], vec_large[k+1][1], 'small', _chr+"_"+str(i)))
                i += 1
            
            vec.append( (_chr, vec_large[k+1][1], vec_large[k+1][2], vec_large[k+1][3], _chr+"_"+str(i)))
            vec.append( (_chr, vec_large[k+1][2], 1000.0 , 'small', _chr+"_"+str(i+1)))
        
        return vec
    
    def write_data(self):
        k= 1
        fh = open(self.data_file, "w")
        for dp in self.data:
            fh.write("fusion_event_"+str(self.sid)+"_dp_"+str(k)+" "+dp[0]._origin.position._chr.replace('chr','hs')+" "+str(dp[0]._origin.position.pos)+" "+str(dp[0]._origin.position.pos+1)+"\n")
            fh.write("fusion_event_"+str(self.sid)+"_dp_"+str(k)+" "+dp[0]._target.position._chr.replace('chr','hs')+" "+str(dp[0]._target.position.pos)+" "+str(dp[0]._target.position.pos+1)+"\n")
            fh.write("\n")
            k += 1
        
        fh.close()
    
    def run(self):
        #1. remove existing circos.png and circos.svg
        for _file in ["circos.png", "circos.svg"]:
            if os.path.exists(_file):
                os.remove(_file)
        
        #2. 
        p = Popen([os.getenv("CIRCOS_DIR")+"/bin/circos", "-conf","share/circos/circos.conf","-debug_group","summary,timer"], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = p.communicate()
        rc = p.returncode
        
        #3. move circos.png and circos.svg into tmp/fusion_event_$id/circos.png etc
        if not os.path.exists("tmp/circos"):
            os.mkdir("tmp/circos")
        
        for _file in ["circos.png", "circos.svg"]:
            os.rename(_file,"tmp/circos/fusion_event_"+str(self.sid)+"_"+_file)
        
