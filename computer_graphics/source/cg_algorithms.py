#!/usr/bin/env python
# -*- coding:utf-8 -*-

# 本文件只允许依赖math库
import math


def draw_line(p_list, algorithm):
    """绘制线段

    :param p_list: (list of list of int(): [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'，此处的'Naive'仅作为示例，测试时不会出现
    :return: (list of list of int(): [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0, x1, y1 = round(p_list[0][0]), round(p_list[0][1]), round(p_list[1][0]), round(p_list[1][1])
    if x0 > x1:
        x0, y0, x1, y1 = x1, y1, x0, y0
    result = []
    if algorithm == 'Naive':
        if x0 == x1:
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            for x in range(x0, x1 + 1):
                result.append((x, round(y0 + k * (x - x0))))
    elif algorithm == 'DDA':
        dx = x1 - x0
        dy = y1 - y0
        xk, yk = x0, y0
        if dx == 0:
            if y0 > y1:
                y0, y1 = y1, y0
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        elif dy == 0:
            for x in range(x0, x1 + 1):
                result.append((x, y0))
        elif abs(dy) <= abs(dx):
            m = dy / dx
            for x in range(abs(dx)+1):
                result.append((xk, round(yk)))
                xk += 1
                yk += m
        else:
            s = 1 if y0 < y1 else -1
            m = dx / dy
            for y in range(abs(dy)+1):
                result.append((round(xk), yk))
                xk += s * m
                yk += s
    elif algorithm == 'Bresenham':
        s = 1 if y0 < y1 else -1
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        xk, yk = x0, y0
        if dy < dx:
            p = 2*dy - dx
            for x in range(dx + 1):
                result.append((xk, yk))
                if p >= 0:
                    yk += s
                    p -= 2*dx
                xk += 1
                p += 2*dy
        else:
            p = 2*dx - dy
            for y in range(dy + 1):
                result.append((xk, yk))
                if p >= 0:
                    xk += 1
                    p -= 2*dy
                yk += s
                p += 2*dx

    return result


def draw_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    for i in range(len(p_list)):
        line = draw_line([p_list[i - 1], p_list[i]], algorithm)
        result += line
    return result

def draw_part_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    for i in range(1, len(p_list)):
        line = draw_line([p_list[i - 1], p_list[i]], algorithm)
        result += line
    return result

def draw_ellipse(p_list):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    quarter = []
    xc, yc = round((x0 + x1) / 2), round((y0+y1)/2)
    rx, ry = round((abs(x1-x0))/2), round((abs(y1-y0))/2)
    p = float(ry**2+rx**2*(0.25-ry))
    xk, yk = 0, ry
    quarter.append((xk, yk))
    while(ry**2*xk<rx**2*yk):
        if(p<0):
            p+=float(ry**2*(2*xk+3))
        else:
            p+=float(ry**2*(2*xk+3)-rx**2*(2*yk-2))
            yk-=1
        xk+=1
        quarter.append((xk, yk))
    p=float((ry*(xk+0.5))**2+(rx*(yk-1))**2-(rx*ry)**2)
    while(yk>0):
        if(p<0):
            p+=float(ry**2*(2*xk+2)+rx**2*(-2*yk+3))
            xk+=1
        else:
            p+=float(rx**2*(-2*yk+3))
        yk-=1
        quarter.append((xk, yk))
    result = []
    for p in quarter:
        result.append((xc+p[0], yc+p[1]))
        result.append((xc+p[0], yc-p[1]))
        result.append((xc-p[0], yc+p[1]))
        result.append((xc-p[0], yc-p[1]))
    return result


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    du = 0.001
    result = []
    n = len(p_list)
    if algorithm == 'Bezier':
        result.append(p_list[0])
        u = du
        while u < 1:
            p = p_list.copy()
            for i in range(n-1):
                p_new = []
                for k in range(len(p) - 1):
                    p_new.append([(1 - u) * p[k][0] + u * p[k+1][0], (1 - u) * p[k][1] + u * p[k+1][1]])
                p = p_new.copy()
            x, y = round(p[0][0]), round(p[0][1])
            result.append([x, y])
            u += du
        result.append(p_list[-1])
    elif algorithm == 'B-spline':
        u = 0
        while u <= 1:
            coef = [-u ** 3 + 3 * u ** 2 - 3 * u + 1, 3 * u ** 3 - 6 * u ** 2 + 4,
                    -3 * u ** 3 + 3 * u ** 2 + 3 * u + 1, u ** 3]
            for i in range(n - 3):
                x, y = 0.0, 0.0
                for j in range(4):
                    x += coef[j] * p_list[i+j][0]
                    y += coef[j] * p_list[i+j][1]
                result.append([round(x/6), round(y/6)])
            u += du
    return result


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    return [[p[0] + dx, p[1] + dy] for p in p_list]


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    theta = r * math.pi / 180
    for p in p_list:
        tx = round((p[0] - x) * math.cos(theta) - (p[1] - y) * math.sin(theta) + x)
        ty = round((p[0] - x) * math.sin(theta) + (p[1] - y) * math.cos(theta) + y)
        p[0], p[1] = tx, ty
    return p_list


def scale(p_list, x, y, s):
    """缩放变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    return [[round((p[0] - x) * s + x), round((p[1] - y) * s + y)] for p in p_list]


def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """线段裁剪

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :param algorithm: (string) 使用的裁剪算法，包括'Cohen-Sutherland'和'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """
    result = []
    if x_min == x_max or y_min == y_max:
        result = [[0, 0], [0, 0]]
        return result
    if x_min > x_max:
        x_min, x_max = x_max, x_min
    if y_min > y_max:
        y_min, y_max = y_max, y_min
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    if algorithm == 'Cohen-Sutherland':
        while True:
            code0, code1 = 0, 0
            if x0 < x_min:
                code0 += 1
            elif x0 > x_max:
                code0 += 2
            if y0 < y_min:
                code0 += 4
            elif y0 > y_max:
                code0 += 8
            if x1 < x_min:
                code1 += 1
            elif x1 > x_max:
                code1 += 2
            if y1 < y_min:
                code1 += 4
            elif y1 > y_max:
                code1 += 8
            if (code0 | code1) == 0:
                result = [[x0, y0], [x1, y1]]
                break
            elif (code0 & code1) != 0:
                result.append([0, 0])
                result.append([0, 0])
                break
            else:
                if code0 == 0:
                    x0, x1 = x1, x0
                    y0, y1 = y1, y0
                    code0, code1 = code1, code0
                if code0 & 1:
                    y0 = round(y0 + ((x_min - x0) * (y0 - y1) / (x0 - x1)))
                    x0 = x_min
                if code0 & 2:
                    y0 = round(y0 + ((x_max - x0) * (y0 - y1) / (x0 - x1)))
                    x0 = x_max
                if code0 & 4:
                    x0 = round(x0 + ((y_min - y0) * (x0 - x1) / (y0 - y1)))
                    y0 = y_min
                if code0 & 8:
                    x0 = round(x0 + ((y_max - y0) * (x0 - x1) / (y0 - y1)))
                    y0 = y_max

    elif algorithm == 'Liang-Barsky':
        q = [x0 - x1, x1 - x0, y0 - y1, y1 - y0]
        d = [x0 - x_min, x_max - x0, y0 - y_min, y_max - y0]
        t0, t1 = 0, 1
        for i in range(4):
            if q[i] < 0:
                t0 = max(t0, d[i] / q[i])
            elif q[i] > 0:
                t1 = min(t1, d[i] / q[i])
            elif q[i] == 0 and d[i] < 0:
                result = [[0, 0], [0, 0]]
                return result
            if t0 > t1:
                result = [[0, 0], [0, 0]]
                return result
        result = [[round(x0 + t0 * (x1 - x0)), round(y0 + t0 * (y1 - y0))],
                  [round(x0 + t1 * (x1 - x0)), round(y0 + t1 * (y1 - y0))]]
    return result