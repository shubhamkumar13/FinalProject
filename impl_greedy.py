from math import sqrt, ceil, floor
import numpy as np
from scipy.spatial import ConvexHull
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row
from bokeh.models import ColumnDataSource, DataRange1d, Plot, LinearAxis, Grid
from bokeh.models.glyphs import Circle, Segment
from random import randint
from matplotlib import collections as mc
import collections as cl
from sympy.geometry import Segment as Line
from sympy.geometry import Point as Point
from sympy.geometry import Circle as Circle
from sympy.geometry import Line as SympyLine
import pprint
import random
from itertools import zip_longest
from sympy import Symbol
import sys

# orig_stdout = sys.stdout
# f = open('/home/shubham/lsdc/output.txt', 'a')
# sys.stdout = f

# class Point(object):
#     def __init__(self, x, y):
#         self.x = x
#         self.y = y
#         self.point = Point(x, y)

#     def __repr__(self):
#         return "( %s, %s)" % (self.x, self.y)


# class Line(object):

#     def __init__(self, p1, p2):
#         self.point1 = p1
#         self.point2 = p2
#         self.line = Segment(p1, p2)

#     def Length(self):
#         l = self.line.length
#         return l

#     def Slope(self):
#         l = SympyLine(self.point1, self.point2)
#         m = l.slope
#         return m

#     def Y_intercept(self):
#         m = self.Slope()
#         return (self.point1.y - m*self.point1.x)

#     def X_intercept(self):
#         return self.point1.x

#     def __repr__(self):
#         m = self.Slope()
#         c = self.Y_intercept()
#         if m == float("inf"):
#             c = self.X_intercept()
#             return "( x = %s )\n" % (c)
#         elif m is 0:
#             return "( y = %s )\n" % (c)
#         else:
#             return "( y = %sx + %s )\n" % (m, c)


# class Circle(object):
#     """docstring for Circle"""

#     def __init__(self, p0, radius):
#         self.point_origin = p0
#         self.radius = radius
#         self.circle = Circle(p0, radius)

#     def __repr__(self):
#         return "( x - %s )^2 + ( y - %s )^2 = %s^2\n" % (self.point_origin.x, self.point_origin.y, radius)


class lineCollection(object):
    '''
        return a dictionary of lines
    '''
    dict_lines = cl.OrderedDict()
    dict_lines_length = cl.OrderedDict()

    def __init__(self, lines):
        self.lines = lines
        self.dict_lines = self.DictLines()
        self.dict_lines_length = self.DictLinesLength()

    def NumLines(self):
        c = len(self.lines)
        return c

    def DictLines(self):
        c = self.NumLines()
        key = [i for i in range(c)]
        value = [line for line in self.lines]
        dict_lines = cl.OrderedDict(zip(key, value))
        # self.dict_lines = dict_lines
        return dict_lines

    def DictLinesLength(self):
        c = self.NumLines()
        key = [i for i in range(c)]
        lengths = [line.length.evalf() for line in self.lines]
        dict_lines_length = cl.OrderedDict(
            zip(key, lengths))
        # self.dict_lines_length = dict_lines_length
        return dict_lines_length

    def __repr__(self):
        disp = list()
        for key, line in self.dict_lines.items():
            disp.append(line.__repr__())

        strval = '\n'.join(disp)
        return strval


class circleCollection(object):
    '''
        return a dictionary of circles
    '''
    dict_circle = cl.OrderedDict()

    def __init__(self, circles):
        self.circles = circles
        self.dict_circles = self.DictCircles()

    def NumCircles(self):
        c = len(self.circles)
        return c

    def DictCircles(self):
        c = self.NumCircles()
        key = [i for i in range(c)]
        value = [circle for circle in self.circles]
        dict_circles = cl.OrderedDict(zip(key, value))
        return dict_circles

    def __repr__(self):
        strval = [val.__repr__() for key, val in self.dict_circles.items()]
        strval = '\n'.join(strval)
        return strval


class Intersection(object):
    '''
        all the collections of lines and circles
        are present inside this
        and all collections of intersections
        are calculated
    '''
    line_and_circle = cl.OrderedDict()
    circle_and_circle = cl.OrderedDict()
    common_dict = cl.OrderedDict()
    line_circle_length = cl.OrderedDict()
    circle_line_length = cl.OrderedDict()
    circle_and_line = cl.OrderedDict()
    delta = cl.OrderedDict()
    new_delta = cl.OrderedDict()
    # circle_line_set = set()
    # circle_circle_set = set()
    # line_circle_set = set()
    circle_pairs = list()

    def __init__(self, dict_lines, dict_circles):
        print("reached init intersection")
        self.dict_lines = dict_lines
        self.dict_circles = dict_circles
        self.circle_and_line = self.CircleLine()
        self.circle_and_circle = self.CircleCircle()
        self.circle_line_length = self.CircleLineSum()
        self.circle_pairs = self.CirclePairList()
        self.common_dict = self.CirclePairLineDict()
        # self.line_circle_length = self.LineCircleSum()
        # self.line_and_circle = self.LineCircle()
        # self.circle_line_set = set(self.circle_and_line)
        # self.circle_circle_set = set(self.circle_and_circle)
        # self.line_circle_set = set(self.line_and_circle)
        print("Bye Init")

    # def Intersects(self, circle1, line=None, circle2=None):
    #     if line is None:
    #         points = circle1.intersection(circle2)
    #         if len(points) > 0:
    #             return True
    #         else:
    #             return False
    #     else:
    #         points = line.intersection(circle1)
    #         if len(points) > 0:
    #             return True
    #         else:
    #             return False

    def TwoCirclesOneLine(self, circle1index, circle2index, line):
        circle_intersection_list = self.circle_and_circle[circle1index]

        if circle2index in circle_intersection_list:

            line_interesection_list = self.circle_and_line[circle2index]

            if line in line_interesection_list:

                line_seg = self.dict_lines[line]
                circle1_val = self.dict_circles[circle1index]
                circle2_val = self.dict_circles[circle2index]
                
                v_unit = dict()
                v_unit = {
                    "p1 out c1" : not circle1_val.encloses_point(line_seg.points[0]),
                    "p2 out c1" : not circle1_val.encloses_point(line_seg.points[1]),
                    "p1 out c2" : not circle2_val.encloses_point(line_seg.points[0]),
                    "p2 out c2" : not circle2_val.encloses_point(line_seg.points[1]),
                    "p1 in c1" : circle1_val.encloses_point(line_seg.points[0]),
                    "p2 in c1" : circle1_val.encloses_point(line_seg.points[1]),
                    "p1 in c2" : circle2_val.encloses_point(line_seg.points[0]),
                    "p2 in c2" : circle2_val.encloses_point(line_seg.points[1])
                }

                # v_binary = {
                #     "p1 and p2 out c1": v_unit["p1 out c1"] and v_unit["p2 out c1"],
                #     "p1 and p2 out c2": v_unit["p1 out c2"] and v_unit["p2 out c2"],
                # }

                # v_ternary = {
                #     "p1 and p2 out c1 not c2": v_binary["p1 and p2 out c1"] and (not(v_unit["p1 out c2"]) or not(v_unit["p2 out c2"])),
                #     "p1 and p2 out c2 not c1": v_binary["p1 and p2 out c2"] and (not(v_unit["p1 out c1"]) or not(v_unit["p2 out c1"])),
                #     "p1 and p2 out c1 and c2": v_binary["p1 and p2 out c2"] and v_binary["p1 and p2 out c1"],
                # }

                case = {
                    "A": (v_unit["p1 out c1"] and v_unit["p2 out c1"] and v_unit["p1 out c2"] and v_unit["p2 out c2"]),

                    "B": ((v_unit["p1 out c1"] and v_unit["p2 out c1"])  
                        and ((v_unit["p2 in c2"] and v_unit["p1 out c2"]) or (v_unit["p2 in c1"] and v_unit["p2 out c1"]))),
                        
                        # or
                        
                        # ((v_unit["p1 out c2"] and v_unit["p2 out c2"])  
                        # and ((v_unit["p1 in c1"] and v_unit["p2 out c1"]) or (v_unit["p2 in c1"] and v_unit["p1 out c1"]))),

                    "C": (v_unit["p1 out c1"] and v_unit["p2 in c1"] and v_unit["p1 out c2"] and v_unit["p2 in c2"])
                        or (v_unit["p1 in c1"] and v_unit["p2 out c1"] and v_unit["p1 in c2"] and v_unit["p2 out c2"]),
                        
                        # or
                        
                        # (((v_unit["p1 in c1"] and v_unit["p2 out c1"]) and (v_unit["p1 in c2"] and v_unit["p2 out c2"]))
                        # or ((v_unit["p2 in c1"] and v_unit["p1 out c1"]) and (v_unit["p2 in c2"] and v_unit["p1 out c2"]))),
                    
                    "D": ((v_unit["p1 in c1"] and v_unit["p2 in c1"]) and (v_unit["p2 in c2"] and v_unit["p1 out c2"]))
                        or ((v_unit["p1 in c1"] and v_unit["p2 out c1"]) and (v_unit["p1 in c2"] and v_unit["p2 in c2"])),

                    #     or

                    #    ((v_unit["p1 in c2"] and v_unit["p2 in c2"]) and (v_unit["p1 in c1"] and v_unit["p2 out c1"]))
                    #     or ((v_unit["p2 in c2"] and v_unit["p1 in c2"]) and (v_unit["p1 in c1"] and v_unit["p2 out c1"]))),
                    
                    "E": (v_unit["p1 in c1"] and v_unit["p2 out c1"] and v_unit["p2 in c2"] and v_unit["p1 out c2"])
                        # or (v_unit["p2 in c1"] and v_unit["p1 out c1"] and v_unit["p1 in c1"] and v_unit["p2 out c1"])
                    
                }

                intersects = case["A"] or case["B"] or case["C"] or case["D"] or case["E"]
                points1 = line_seg.intersection(circle1_val)
                points2 = line_seg.intersection(circle2_val)

                if intersects:
                    if case["A"]:
                        p1 = points1[1]
                        p2 = points2[0]
                        if circle1_val.encloses_point(p1) and circle2_val.encloses_point(p2):
                            temp = p1.distance(p2).evalf()
                            return temp
                        else:
                            return None
                    elif case["B"]:
                        if len(points1) == 1:
                            p1 = points1[0]
                            p2 = points2[0]
                            temp = p1.distance(p2).evalf()
                            return temp
                        elif len(points2) == 1:
                            p1 = points1[1]
                            p2 = points2[0]
                            temp = p1.distance(p2).evalf()
                            return temp
                    elif case["C"]:
                        if v_unit["p2 in c1"] and v_unit["p2 in c2"]:
                            p2 = line_seg.points[1]
                            p1 = points2[0]
                            temp = p2.distance(p1).evalf()
                            return temp
                        elif v_unit["p1 in c1"] and v_unit["p1 in c2"]:
                            p1 = line_seg.points[0]
                            p2 = points1[0]
                            temp = p2.distance(p1).evalf()
                            return temp
                    elif case["D"]:
                        if len(points1) == 0:
                            p1 = points2[0]
                            p2 = line_seg.points[1]
                            temp = p2.distance(p1).evalf()
                            return temp
                        elif len(points2) == 0:
                            p1 = line_seg.points[0]
                            p2 = points1[0]
                            temp = p2.distance(p1).evalf()
                            return temp
                    elif case["E"]:
                        p1 = points1[0]
                        p2 = points2[0]
                        temp = p2.distance(p1).evalf()
                        return temp
                else:
                    return None               


                # index = line_interesection_list.index(line)
                # temp_list = self.circle_line_length[circle2index]
                # delta_value2 = temp_list[index]

                # line_list = self.circle_and_line[circle1index]
                # index = line_list.index(line)
                # temp_list = self.circle_line_length[circle1index]
                # delta_value1 = temp_list[index]

                # extra_length = delta_value1 + delta_value2
                # original_length =
                # excess = extra_length - original_length

                # print(delta_value1, delta_value2, original_length)
                # print(excess)
                # return excess

    def CirclePairList(self):
        pair_list = [[key,v] for key, value in self.circle_and_circle.items() for v in value]
        ctr = cl.Counter(frozenset(x) for x in pair_list)
        # boolean = [ctr[frozenset(x)] > 0 for x in pair_list]

        pair_list = [list(k) for k,v in ctr.items()]
        print("Pairs")
        pprint.pprint(pair_list)
        print("\n")
        return pair_list
    
    def overlap(self, min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))
    
    def checkval(self, val):
        if val < 1:
            return 0
        else:
            return val

    def GetCommonVal(self, line_tuple1, line_tuple2, circle_pair):
        val = float()
        if line_tuple1[0] == line_tuple2[0]:
            p1 = line_tuple1[1]
            p2 = line_tuple2[1]
            line = self.dict_lines[line_tuple1[0]]
            t1, t2 = line.points
            dist = line.length
            if len(p1) == 2 and len(p2) == 2:
                a1, a2 = p1[0], p1[1]
                b1, b2 = p2[0], p2[1]
                
                norm_a1 = a1.distance(t1) / dist
                norm_a2 = a2.distance(t1) / dist
                norm_b1 = b1.distance(t1) / dist
                norm_b2 = b2.distance(t1) / dist

                norm = [norm_a1, norm_a2, norm_b1, norm_b2]
                norm.sort()
                val = self.overlap(norm[0], norm[2], norm[1], norm[3])
                val = float(val*dist)

                # # if a1.distance(b2) >= a2.distance(b1):
                # #     val = float(a2.distance(b1))
                # # if a1.distance(b2) <= a2.distance(b1):
                # #     val = float(a1.distance(b2))
                # d1 = line_tuple1[2]
                # d2 = line_tuple2[2]
                # if a1.distance(a2) < a1.distance(b2):
                #     val = float(b1.distance(a2))
                # elif a1.distance(a2) > a1.distance(b2):
                #     val = float(b1.distance(b2))
                # elif b1.distance(b2) > a1.distance(b2):
                #     val = float(a1.distance(a2))
                    
            elif len(p1) == 2 and len(p2) == 1:
                a1, a2 = p1[0], p1[1]
                b1 = p2[0]
                norm_a1 = a1.distance(t1) / dist
                norm_a2 = a2.distance(t1) / dist
                norm_b1 = b1.distance(t1) / dist
                norm_b2 = 1.0
                    

                norm = [norm_a1, norm_a2, norm_b1, norm_b2]
                norm.sort()
                val = self.overlap(norm[0], norm[2], norm[1], norm[3])
                val = float(val*dist)
                # if a1.distance(a2) > a1.distance(b1) or not(a1.distance(a2) < a2.distance(b1)):
                #     val = float(a2.distance(b1))
                # elif a1.distance(a2) < b1.distance(a2) or a1.distance(a2) < a2.distance(b1):
                
                #     val = float(a1.distance(a2)) 
            elif len(p1) == 1 and len(p2) == 2:
                a2 = p1[0]
                b1,b2 = p2[0], p2[0]
                
                norm_a1 = 0
                norm_a2 = a2.distance(t1) / dist
                norm_b1 = b1.distance(t1) / dist
                norm_b2 = b2.distance(t1) / dist
            
                
                norm = [norm_a1, norm_a2, norm_b1, norm_b2]
                norm.sort()
                val = self.overlap(norm[0], norm[2], norm[1], norm[3])
                val = float(val*dist)

                # if b1.distance(b2) > b1.distance(a1):
                #     val = float(a1.distance(b1))
                # elif b2.distance(b1) < a1.distance(b2) or b2.distance(b1) < a1.distance(b1):
                #     val = float(a1.distance(b2))
            elif len(p1) == 1 and len(p2) == 1:
                a1, b1 = p1[0], p2[0]
                c1, c2 = self.dict_circles[circle_pair[0]], self.dict_circles[circle_pair[1]]

                if c1.encloses_point(t1) and c2.encloses_point(t1) and not(c1.encloses_point(t2)) and not(c2.encloses_point(t1)):
                    val = min(a1.distance(t1).evalf(), b1.distance(t1).evalf())
                    val = float(val)
                if c1.encloses_point(t2) and c2.encloses_point(t2) and not(c1.encloses_point(t1)) and not(c2.encloses_point(t1)):
                    val = min(a1.distance(t2).evalf(), b1.distance(t2).evalf())
                    val = float(val)
                   
                # if c2.encloses_point(a1) and c2.encloses_point(t1) and not(c2.encloses_point(t2)):
                #     val = float(t1.distance(a1))
                # elif c1.encloses_point(b1) and c1.encloses_point(t2) and not(c1.encloses_point(t1)):
                #     val = float(t2.distance(b1))
                # elif (c1.encloses_point(b1) and c1.encloses_point(t1)) and (c2.encloses_point(a1) and c2.encloses_point(t2)):
                #     val = float(a1.distance(b1))

        else:
            val = None

        val = self.checkval(val)
        return val
        
    def CirclePairLineDict(self):
        temp_dict = dict()
        for pair in self.circle_pairs:
            list1 = self.circle_and_line[pair[0]]
            list2 = self.circle_and_line[pair[1]]
            set1 = set(list1)
            set2 = set(list2)
            s = list(set1.intersection(set2))
            print(pair)
            for element in s:
                print(element)
                print(list1)
                print(list2)
                index1 = list1.index(element)
                index2 = list2.index(element)
                len_list1 = self.circle_line_length[pair[0]]
                len_list2 = self.circle_line_length[pair[1]]
                value = self.GetCommonVal(len_list1[index1], len_list2[index2], pair)
                pair = tuple(pair)
                if value is not None and value != 0:
                    temp_dict[(pair, element)] = value
        
        return temp_dict

            

        

    

    # def GetCommonDict(self, circle):
    #     line_interesection_list=self.circle_and_line[circle]
    #     circle_intersection_list=self.circle_and_circle[circle]
    #     for circle2 in circle_intersection_list:
    #         for line in line_interesection_list:
    #             value_common=self.TwoCirclesOneLine_temp(circle, circle2, line)
    #             if value_common is not None:
    #                 if float(value_common) > 0.0 and float(value_common) < 1.0:
    #                     self.common_dict[frozenset([circle, circle2]), line] = floor(value_common) 
    #                 else:
    #                     self.common_dict[frozenset([circle, circle2]), line]= value_common  


    # def GenerateDisks(self):
    #     print("Creating new disks")

    #     val = [value for key,value in self.delta.items()]

    #     for key,value in self.delta.items():
    #         self.RemoveDelta(key)





    # The following three functions create
    # indexed values of intersections
    # 1-- Circle intersects Line -------> Circle(index) -> [Collection of Line(index) who intersect Circle(index)]
    # 2-- Circle intersects Circle -----> Circle(index) -> [Collection of Circle(index) who intersect Circle(index)]
    # 3-- Line intersects Circle -------> Line(index) -> [Collection of Circle(index) who intersect Line(index)]

    def CircleLine(self):
        print("reached circle and line dict creation")
        temp_list=list()
        temp_dict=dict()
        for key_circle, circle in self.dict_circles.items():
            for key_line, line in self.dict_lines.items():
                points=circle.intersection(line)
                if len(points) > 0:
                    temp_list.append(key_line)
            temp_dict[key_circle]=temp_list
            temp_list=[]
        temp_dict=cl.OrderedDict(temp_dict)
        print("bbye circle and line dict creation")
        return temp_dict

    def CircleCircle(self):
        print("reached circle and circle dict creation")
        temp_list=list()
        temp_dict=dict()
        for key_circle_main, circle_main in self.dict_circles.items():
            for key_circle, circle in self.dict_circles.items():
                if key_circle != key_circle_main:
                    points=circle_main.intersection(circle)
                    if len(points) == 2 and points[0].distance(points[1]) > 1:
                        temp_list.append(key_circle)
            temp_dict[key_circle_main]=temp_list
            temp_list=list()
        temp_dict=cl.OrderedDict(temp_dict)
        print("bbye circle and circle dict creation")
        return temp_dict

    # def LineCircle(self):
    #     print("reached line and circle dict creation")
    #     temp_list = list()
    #     temp_dict = dict()
    #     for key_line, line in self.dict_lines.items():
    #         for key_circle, circle in self.dict_circles.items():
    #             points = circle.intersection(line)
    #             if len(points) > 0:
    #                 temp_list.append(key_circle)
    #         temp_dict[key_line] = temp_list
    #         temp_list = []
    #     temp_dict = cl.OrderedDict(temp_dict)
    #     print("end line and circle dict")
    #     return temp_dict
    # The Following function gives out the length of intersection
    # Only when They are SYMPY objects
    # Amount line present in a circle

    # def LengthIntersection(self, line, circle):
    #     p1, p2 = line.points
    #     dist_p1, dist_p2 = p1.distance(
    #         circle.center), p2.distance(circle.center)

    #     if dist_p1 < radius:
    #         intersection_point = circle.intersection(line)
    #         dist_segment = intersection_point[0].distance(p1)
    #         return dist_segment

    #     elif dist_p2 < radius:
    #         intersection_point = circle.intersection(line)
    #         dist_segment = intersection_point[0].distance(p2)
    #         return dist_segment

    #     else:
    #         intersection_point = circle.intersection(line)
    #         dist_segment = intersection_point[0].distance(
    #             intersection_point[1])
    #         return dist_segment


    # get individual circle line intersection
    # then calculate individual length of each
    # circle and the respective lengths
    # make a dictionary of
    # circle_index -> [list of linesum present inside]
    # circle_index -> [total sum of all lines]

    def GetIndividualCircleLength(self, index):
        length_segment=list()
        line_index=self.circle_and_line[index]
        circle=self.dict_circles[index]
        for line in line_index:

            line_seg=self.dict_lines[line]

            temp=line_seg.intersection(circle)
            if circle.encloses_point(line_seg.p1):
                temp_point=temp[0]
                temp_val=line_seg.p1.distance(temp[0]).evalf()
                temp[0] = Point(float(temp_point.x), float(temp_point.y))
            elif circle.encloses_point(line_seg.p2):
                temp_point=temp[0]
                temp_val=line_seg.p2.distance(temp_point).evalf()
                temp[0] = Point(float(temp_point.x), float(temp_point.y))
            elif circle.encloses_point(line_seg.p1) and circle.encloses_point(line_seg.p1):
                temp_point1=line_seg.p1
                temp_point2=line_seg.p2
                temp_val=temp_point2.distance(temp_point1).evalf()
                temp[0] = Point(float(temp_point1.x), float(temp_point1.y))
                temp[1] = Point(float(temp_point2.x), float(temp_point2.y))
            elif not(circle.encloses_point(line_seg.p1) and circle.encloses_point(line_seg.p2)):
                if len(temp) == 2:
                    temp_point1=temp[0]
                    temp_point2=temp[1]
                    temp_val=temp_point2.distance(temp_point1).evalf()
                    temp[0] = Point(float(temp_point1.x), float(temp_point1.y))
                    temp[1] = Point(float(temp_point2.x), float(temp_point2.y))
                else:
                    continue

            temp_val = float(temp_val)
            if temp_val < 1:
                temp_val = 0

            length_segment.append((line,temp,temp_val))

        return length_segment

    def ListSum(self, list_val):
        sum=0
        for (line,points,value) in list_val:
            sum=sum + value
        return sum

    def CircleLineSum(self):
        delta_dict=cl.OrderedDict()
        temp_dict=cl.OrderedDict()
        delta_sum=0
        for key, value in self.dict_circles.items():
            delta_dict[key]=self.GetIndividualCircleLength(key)

        for key, value in delta_dict.items():
            delta_sum=self.ListSum(value)
            temp_dict[key]= delta_sum

        temp_dict=sorted(temp_dict.items(), key=lambda t: t[1], reverse=True)
        self.delta=temp_dict
        self.delta=cl.OrderedDict(self.delta)
        return delta_dict


    # get individual line length
    # and make 2 dictionaries
    # 1---line_index -> circle_intersection_index
    # 2---line_index -> circle_intersection_value


    # def GetIndividualLineLength(self, line_index, circle_index):
    #     line = self.dict_lines[line_index]
    #     circle = self.dict_circles[circle_index]

    #     points = line.intersection(circle)
    #     assert (len(points) > 0), "no points intersect, error in line's sum"
    #     if circle.encloses_point(line.p1):
    #         temp_val = line.p1.distance(points[0])
    #     elif circle.encloses_point(line.p2):
    #         temp_val = line.p2.distance(points[0])
    #     else:
    #         temp_val = points[1].distance(poins[0])

    #     return temp_val

    # def FindAllCircleSum(self, line_index, circle_list):
    #     line_delta = list()
    #     for circle_index in circle_list:
    #         temp = self.GetIndividualLineLength(line_index, circle_index)
    #         line_delta.append(temp)
    #     return line_delta

    # def LineCircleSum(self):
    #     temp_dict = cl.OrderedDict()

    #     for key, value in self.line_and_circle.items():
    #         temp_dict[key] = self.FindAllCircleSum(key, value)

    #     self.line_circle_length = temp_dict

    def GetNewDelta(self):
        temp_delta = cl.OrderedDict(self.delta)
        delta_key = [k for k,v in temp_delta.items()]
        key_list = [k for k,v in temp_delta.items()]
        common_dict = self.common_dict
        
        for k, v in common_dict.items():
            c, l = k
            c1, c2 = c
            if delta_key.index(c1) > delta_key.index(c2):
                temp = temp_delta[c1]
                temp = temp - v
                temp = self.checkval(temp)
                temp_delta[c1] = temp
            if delta_key.index(c1) < delta_key.index(c2):
                temp = temp_delta[c2]
                temp = temp - v
                temp = self.checkval(temp)
                temp_delta[c2] = temp

        
        for key in key_list:
            if temp_delta[key] >= 1:
                self.new_delta[key] = temp_delta[key]

        # for key in key_list:
        #     self.new_delta[key] = temp_delta.pop(key, None)
        #     for k, v in common_dict.items():
        #         c, l = k
        #         if key == c[0] and c[1] in temp_delta:
        #             temp = float(temp_delta[c[1]])
        #             temp = abs(temp - v)
        #             temp = self.checkval(temp)
        #             temp_delta[c[1]] = temp


    

def disk_removal(circle_obj, line_obj):
    pass

def disk_cover_algorithm(circle_obj, line_obj):
    dict_lines=line_obj.dict_lines
    dict_circles=circle_obj.dict_circles
    intersection_obj=Intersection(dict_lines, dict_circles)
    # old_intersection_obj = Intersection(dict_lines, dict_circles)

    line_intersection_dict=cl.OrderedDict()
    line_intersection_dict=intersection_obj.circle_and_line
    circle_intersection_dict=intersection_obj.circle_and_circle

    print("\n\n")
    print("Indexed Lines Generated")
    pprint.pprint(dict_lines)
    print("\n\n")    
    print("Indexed Circles Generated")
    pprint.pprint(dict_circles)
    print("\n\n")
    print("Circle and Line Intersection Index")
    pprint.pprint(line_intersection_dict)
    print("\n\n")
    print("Circle and Circle Interesection Index")    
    pprint.pprint(circle_intersection_dict)
    print("\n\n")    

    dict_delta=cl.OrderedDict()
    delta=cl.OrderedDict()
    dict_delta=intersection_obj.circle_line_length
    # old_delta = old_intersection_obj.circle_line_length
    delta=intersection_obj.delta
    print("\n\n")
    print("Circle and Length(Line) Interesection Index")    
    pprint.pprint(dict_delta)
    print("\n\n")

    print("\n\n")
    print("Circles and Sum(Length(Line) Intersection Index")        
    pprint.pprint(delta)
    print("\n\n")    

    # print("Common Dictionary")
    # pprint.pprint(common_dict)
    # print("\n\n")
    
    common_dict = intersection_obj.common_dict
    print("Common Dict \n")
    pprint.pprint(common_dict)
    print("\n")
    intersection_obj.GetNewDelta()
    new_delta = intersection_obj.new_delta
    new_delta = cl.OrderedDict(new_delta)

    print("New reduced delta")
    pprint.pprint(new_delta)
    print("\n\n")

    key_old = [k for k,v in delta.items()]
    key_new = [k for k,v in new_delta.items() if v > 0]
    non_redundant_circles = set(key_old).intersection(key_new)
    non_redundant_circles = list(non_redundant_circles)
    
    print("Non Redundant Disks")
    pprint.pprint(non_redundant_circles)
    print("\n\n")

    return non_redundant_circles
    



def grouper(iterable, n, fillvalue=None):
    # Collect data into fixed-length chunks or blocks
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args=[iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def Random_Lines(points):
    '''
        input points in pairs which form lines
        since points are genrated randomly
        need to take point pairs
        that would serve as input lines
    '''
    Lines=list()
    random.shuffle(points)
    # pprint.pprint(points)
    for p1, p2 in grouper(points, 2):
        line=Line(p1, p2)
        Lines.append(line)
    return Lines

def FindEqPoints(line, num, radius):
    temp_list=list()

    p=line.points
    init_point=p[0]
    exit_point=p[1]
    # m = -(coeff[0]/coeff[1])
    # c = -(coeff[2]/ coeff[1])
    dia=radius*2
    d=radius*2
    length=line.length.evalf()

    const1=(exit_point.x - init_point.x)/length
    const2=(exit_point.y - init_point.y)/length

    x=init_point.x
    y=init_point.y

    num=int(num)
    for i in range(num):

        x=x + const1*d
        y=y + const2*d

        temp_list.append(Point(x, y))
    d=d + dia
    return temp_list

def convertToCircles(temp_list, line, radius):
    p=line.points
    init_point=p[0]
    circles=list()

    for point in temp_list:
        p0=init_point.midpoint(point)
        circles.append(Circle(p0, radius))
        init_point=point

    return circles

def Planned_Circles(line_obj, radius=5):
    '''
        input points not pairs of points
        this means generate points
        then the unit radii would take care
    '''
    circles=cl.OrderedDict()
    num_circles=list()
    lines=cl.OrderedDict()
    length_lines=cl.OrderedDict()


    lines=line_obj.dict_lines

    length_lines=line_obj.dict_lines_length
    # pprint.pprint(length_lines)
    dia=radius * 2

    for key, val in lines.items():
        len_line=length_lines[key]
        num=len_line / dia
        print(len_line)
        num=ceil(num)
        print(num)
        num_circles.append(num)

    for key, line in lines.items():
        c=num_circles[key]
        p_list=FindEqPoints(line, c, radius)
        pprint.pprint(p_list)
        temp_list=convertToCircles(p_list, line, radius)
        circles[key]=temp_list

    circle_list=list()
    for key, val in circles.items():
        circle_list=circle_list + val

    temp_list=list()
    for circle in circle_list:
        p0=circle.center
        r=circle.radius
        temp=Circle(p0, r)
        temp_list.append(temp)

    return temp_list

def Unplanned_Circles(points, radius):
    circles = [Circle(p, radius) for p in points]
    return circles


def display_encompass_points(points):
    '''
        As stated we need to surround
        the total lines generated
        so we need a 2d plane of sorts
        this means we need to apply
        Convex Hull
        And we can get back the vertices
        which are part of the convex hull
    '''
    hull=ConvexHull(points)

    import matplotlib.pyplot as plt
    plt.plot(points[:, 0], points[:, 1], '+')
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

    plt.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'r--', lw=2)
    plt.plot(points[hull.vertices[0], 0], points[hull.vertices[0], 1])
    plt.show()


def display_lines(lines):
    '''
        create a display of lines
    '''


def display_circles(circles):
    '''
        create a display of circles
    '''
    pass


def DistPointPoint(p1, p2):
    return p1.distance(p2).evalf()


def square(x):
    return x**2


def DistPointLine(point: Point, line: Line):
    '''
    distance between point and line
    | ax + by + c |     | mx - y + c |
    --------------- = ------------------
    sqrt(a^2 + b^2)      sqrt(m^2 + 1)

    where x and y are points
    and |-| means absolute value
    '''
    dist=line.distance(point).evalf()

    return dist


def PointInCircle(point, circle):
    '''
    if (x - a)^2 + (y - b)^2 - r^2 < 0 :
        true
    else:
        false

    '''
    dist=DistPointPoint(point, circle.point_origin)
    dist=square(dist)
    radius=circle.radius
    radius=square(radius)
    if dist - radius < 0:
        return True
    else:
        return False




def IsIntersection(circle1: Circle, circle2: Circle=None, line: Line=None):
    '''
        create Sympy objects of circles and line
        create a dict for (index of circle) -> [all(index of line)]
    '''

    if line is None and circle2 is not None:
        points=circle1.intersection(circle2)
        if len(points) > 0:
            return True
    elif line is not None and circle2 is None:

        points=circle1.intersection(line)
        if len(points) > 0:
            return True
    else:
        return False




'''
def convert_list(points, index):
    x = points[index]
    x = np.asarray(x)
    print(x)
    return x
'''
'''
def convert_ndarray(x, y):
    print(x)
    print(y)
'''


def convert_points_matrix(points):
    '''
    convert the points object into
    a 2d numpy array
    '''
    temp_list=list()
    matrix=list()
    for p in points:
        temp_list.append(p.x)
        temp_list.append(p.y)
        matrix.append(temp_list)
        temp_list=[]

    matrix=np.array(matrix)
    return matrix


# Random Points generated using poisson's distribution
def lonely(p, X, r):
    m=X.shape[1]
    x0, y0=p
    x=y=np.arange(-r, r)
    x=x + x0
    y=y + y0

    u, v=np.meshgrid(x, y)

    u[u < 0]=0
    u[u >= m]=m-1
    v[v < 0]=0
    v[v >= m]=m-1

    return not np.any(X[u[:], v[:]] > 0)


def Random_Points(m, r, k):
    '''
        this creates random points,
        the random points are returned
        as an np array
    '''
    # m = extent of sample domain
    # r = minimum distance between points
    # k = samples before rejection
    active_list=[]

    # step 0 - initialize n-d background grid
    X=np.ones((m, m))*-1

    # step 1 - select initial sample
    x0, y0=np.random.randint(0, m), np.random.randint(0, m)
    active_list.append((x0, y0))
    X[active_list[0]]=1

    # step 2 - iterate over active list
    while active_list:
        i=np.random.randint(0, len(active_list))
        rad=np.random.rand(k)*r+r
        theta=np.random.rand(k)*2*np.pi

        # get a list of random candidates within [r,2r] from the active point
        candidates=np.round(
            (rad*np.cos(theta)+active_list[i][0], rad*np.sin(theta)+active_list[i][1])).astype(np.int32).T

        # trim the list based on boundaries of the array
        candidates=[(x, y) for x, y in candidates if x >=
                      0 and y >= 0 and x < m and y < m]

        for p in candidates:
            if X[p] < 0 and lonely(p, X, r):
                X[p]=1
                active_list.append(p)
                break
        else:
            del active_list[i]

    '''
    Converts N x 2 array
    to
    A point class
    '''
    X=np.where(X > 0)
    points=list()
    rangeX=X[0]
    rangeY=X[1]

    for x, y in zip(rangeX, rangeY):
        temp_point=Point(x, y)
        points.append(temp_point)
    return points


    

Points_Lines=Random_Points(70, 20, 10)
Points_Circles = Random_Points(70, 20, 10)

n = len(Points_Lines)
Lines=Random_Lines(Points_Lines)
line_obj=lineCollection(Lines)
radius=5
Circles=Planned_Circles(line_obj, radius)
Circles2=Unplanned_Circles(Points_Circles, radius)
Circles.extend(Circles2)
circle_obj=circleCollection(Circles)
a, b=line_obj.NumLines(), circle_obj.NumCircles()

abc = line_obj.dict_lines

x0 = [float(v.p1.x) for k,v in abc.items()]
x1 = [float(v.p2.x) for k,v in abc.items()]
y0 = [float(v.p1.y) for k,v in abc.items()]
y1 = [float(v.p2.y) for k,v in abc.items()]

xyz = circle_obj.dict_circles

x = [float(c.center.x) for k,c in xyz.items()]
y = [float(c.center.y) for k,c in xyz.items()]

p = figure()

p.annulus(x, y, radius,radius+0.1, color = 'red', line_width = 3, line_color = "black")
p.segment(x0, y0, x1, y1, line_width = 3)

non_redundant_disks = disk_cover_algorithm(circle_obj, line_obj)
# disk_removal(circle_obj, line_obj)

xx = [float(c.center.x) for k,c in circle_obj.dict_circles.items() if k in non_redundant_disks]
yy = [float(c.center.y) for k,c in circle_obj.dict_circles.items() if k in non_redundant_disks]

pp = figure()
pp.annulus(xx, yy, radius,radius+0.1, color = 'red', line_width = 3, line_color = "black")
pp.segment(x0, y0, x1, y1, line_width = 3)

output_file("out.html")

show(row(p,pp))

# sys.stdout = orig_stdout
# f.close()