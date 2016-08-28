from ete3.treeview.faces import StaticItemFace
from PyQt4.QtGui import (QGraphicsRectItem, QGraphicsLineItem,
                         QGraphicsPolygonItem, QGraphicsEllipseItem,
                         QPen, QColor, QBrush, QPolygonF, QFont,
                         QPixmap, QFontMetrics, QPainter,
                         QRadialGradient, QGraphicsSimpleTextItem, QGraphicsTextItem,
                         QGraphicsItem)
from PyQt4.QtCore import Qt,  QPointF, QRect, QRectF

class RulerFace(StaticItemFace):
    """
    To draw plots, usually correlated to columns in alignment

    :argument values : a list of values
    :argument None errors : a list of errors associated to each value. elements of the list can contain a list with lower and upper error, if they are different.
    :argument None colors : a list of colors associated to each value
    :argument None header : a title for the plot
    :argument bar kind : kind of plot, one of bar, curve or sticks.
    :argument None fsize : font size for header and labels
    :argument 100 height : height of the plot (excluding labels)
    :argument None hlines : list of y values of horizontal dashed lines to be drawn across plot
    :argument None hlines_col: list of colors associated to each horizontal line
    :argument None col_width : width of a column in the alignment
    :argument red error_col : color of error bars
    """
    def __init__(self, values, errors=None, colors=None, header='',
                 fsize=9, height = 0, hlines=None, kind='bar',
                 hlines_col = None, extras=None, col_width=11,
                 ylim=None, xlabel='', ylabel=''):

        self.col_w = float(col_width)
        self.height = height
#        self.values = [float(v) for v in values]
        self.values = [int(v) for v in values]
        self.width = self.col_w * len (self.values)
        self.errors = errors if errors else []
        self.colors = colors if colors else ['gray'] * len(self.values)
        self.header = header
        self.fsize = fsize
        if ylim:
            self.ylim = tuple((float(y) for y in ylim))
        else:
            self.ylim = (int(min(self.values)-0.5), int(max(self.values)+0.5))
        self.xlabel = xlabel
        self.ylabel = ylabel

        if self.errors:
            if type(self.errors[0]) is list or type(self.errors[0]) is tuple:
                self._up_err = [float(e[1]) for e in self.errors]
                self._dw_err = [float(-e[0]) for e in self.errors]
            else:
                self._up_err = [float(e) for e in self.errors]
                self._dw_err = [float(-e) for e in self.errors]
        if kind == 'bar':
            self.draw_fun = self.draw_bar
        elif kind == 'stick':
            self.draw_fun = self.draw_stick
        elif kind == 'curve':
            self.draw_fun = self.draw_curve
        else:
            raise 'kind %s not yet implemented... ;)'

        self.hlines = [float(h) for h in hlines] if hlines else [1.0]
        self.hlines_col = hlines_col if hlines_col else ['black']*len(self.hlines)

        self.extras = extras if extras else ['']
        if len (self.extras) != len (self.values):
            self.extras = ['']

        super(RulerFace,
              self).__init__(QGraphicsRectItem(-40, 0, self.width+40,
                                                     self.height+50))
        self.item.setPen(QPen(QColor('white')))

    def update_items(self):
        # draw lines
        #for line, col in zip(self.hlines, self.hlines_col):
        #    self.draw_hlines(line, col)
        # draw plot
        #width = self.col_w
        #for i, val in enumerate(self.values):
        #    self.draw_fun(width * i + self.col_w / 2 , val, i)
        # draw error bars
        #if self.errors:
        #    for i in range(len(self.errors)):
        #        self.draw_errors(width * i + self.col_w / 2 , i)
        # draw x axis
        self.draw_x_axis()
        # draw y axis
        #self.draw_y_axis()
        # put header
        #self.write_header()

    def write_header(self):
        text = QGraphicsSimpleTextItem(self.header)
        text.setFont(QFont("Arial", self.fsize))
        text.setParentItem(self.item)
        text.setPos(0, 5)

    def draw_y_axis(self):
        lineItem = QGraphicsLineItem(0, self.coordY(self.ylim[0]),
                                     0, self.coordY(self.ylim[1]),
                                     parent=self.item)
        lineItem.setPen(QPen(QColor('black')))
        lineItem.setZValue(10)
        max_w = 0
        for y in set(self.hlines + list(self.ylim)):
            lineItem = QGraphicsLineItem(0, self.coordY(y),
                                               -5, self.coordY(y),
                                               parent=self.item)
            lineItem.setPen(QPen(QColor('black')))
            lineItem.setZValue(10)
            text = QGraphicsSimpleTextItem(str(y))
            text.setFont(QFont("Arial", self.fsize-2))
            text.setParentItem(self.item)
            tw = text.boundingRect().width()
            max_w = tw if tw > max_w else max_w
            th = text.boundingRect().height()
            # Center text according to masterItem size
            text.setPos(-tw - 5, self.coordY(y)-th/2)
        if self.ylabel:
            text = QGraphicsSimpleTextItem(self.ylabel)
            text.setFont(QFont("Arial", self.fsize-1))
            text.setParentItem(self.item)
            text.rotate(-90)
            tw = text.boundingRect().width()
            th = text.boundingRect().height()
            # Center text according to masterItem size
            text.setPos(-th -5-max_w, tw/2+self.coordY(sum(self.ylim)/2))

    def draw_x_axis(self):
        #lineItem = QGraphicsLineItem(self.col_w/2,
        #                                   self.coordY(self.ylim[0])+2,
        #                                   self.width-self.col_w/2,
        #                                   self.coordY(self.ylim[0])+2,
        #                                   parent=self.item)
        #lineItem.setPen(QPen(QColor('black')))
        #lineItem.setZValue(10)
        #all_vals = list(range(0, len(self.values), 5))
        #if (len(self.values)-1)%5:
        #    all_vals += [len(self.values)-1]
        for x, lab in enumerate(self.values):
#            lineItem = QGraphicsLineItem(0, self.coordY(self.ylim[0])+2,
#                                               0, self.coordY(self.ylim[0])+6,
#                                               parent=self.item)
#            lineItem.setX(x*self.col_w + self.col_w/2)
#            lineItem.setPen(QPen(QColor('black')))
#            lineItem.setZValue(10)
            text = QGraphicsSimpleTextItem(str(lab))
            text.rotate(-90)
            text.setFont(QFont("Arial", self.fsize-2))
            text.setParentItem(self.item)
            tw = text.boundingRect().height()
            # Center text according to masterItem size
            text.setPos(x*self.col_w-tw/2 + self.col_w/2,
                        self.coordY(self.ylim[0]))
  #+6)

    def coordY(self, y):
        """
        return the transformation of Y according to mean value
        (that is last element of lines)
        """
        y_offset = 30
        if self.ylim[1] <= y: return y_offset
        if self.ylim[1] == 0: return self.height + y_offset
        if self.ylim[0] >= y: return self.height + y_offset
        #return self.height - y * self.height / self.ylim[1]
        return self.height + y_offset - (y-self.ylim[0]) / (self.ylim[1]-self.ylim[0]) * self.height

    def draw_hlines (self, line, col):
        lineItem = QGraphicsLineItem(0, self.coordY(line),
                                           self.width, self.coordY(line),
                                           parent=self.item)
        lineItem.setPen(QPen(QColor(col), 1, Qt.DashLine))
        lineItem.setZValue(10)

    def draw_bar(self, x, y, i):
        h = self.coordY(self.ylim[0])#self.height
        coordY = self.coordY
        item = self.item
        # if value stands out of bound
        if y < self.ylim[0]: return
        if y < self.ylim[1]:
            # left line
            lineItem = QGraphicsLineItem(0, h, 0, coordY(y), parent=item)
            lineItem.setX(x-3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))
            # right line
            lineItem = QGraphicsLineItem(0, h, 0, coordY(y), parent=item)
            lineItem.setX(x+3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))
            # top line
            lineItem = QGraphicsLineItem(0, coordY(y), 6, coordY(y), parent=item)
            lineItem.setX(x-3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))
        else:
            # lower left line
            lineItem = QGraphicsLineItem(0, h, 0, coordY(y), parent=item)
            lineItem.setX(x-3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))
            # lower right line
            lineItem = QGraphicsLineItem(0, h, 0, coordY(y), parent=item)
            lineItem.setX(x+3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))
            # upper left line
            lineItem = QGraphicsLineItem(0, coordY(y)-4, 0, coordY(y)-7, parent=item)
            lineItem.setX(x-3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))
            # upper right line
            lineItem = QGraphicsLineItem(0, coordY(y)-4, 0, coordY(y)-7, parent=item)
            lineItem.setX(x+3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))
            # top line
            lineItem = QGraphicsLineItem(0, coordY(y)-7, 6, coordY(y)-7, parent=item)
            lineItem.setX(x-3)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))

    def draw_stick(self, x, y, i):
        lineItem = QGraphicsLineItem(0, self.coordY(self.ylim[0]),
                                           0, self.coordY(y),
                                           parent=self.item)
        lineItem.setX(x)
        lineItem.setPen(QPen(QColor(self.colors[i]),2))

    def draw_errors(self, x, i):
        lower = self.values[i]+self._dw_err[i]
        upper = self.values[i]+self._up_err[i]
        lineItem = QGraphicsLineItem(0, self.coordY(lower), 0,
                                           self.coordY(upper), parent=self.item)
        lineItem.setX(x)
        lineItem.setPen(QPen(QColor('black'),1))

    def draw_curve(self, x, y, i):
        # top line
        lineItem = QGraphicsLineItem(0, self.coordY(y), 4,
                                           self.coordY(y), parent=self.item)
        lineItem.setX(x-2)
        lineItem.setPen(QPen(QColor(self.colors[i]),2))
        if i > 0:
            prev = self.values[i-1] if i>0 else self.values[i]
            lineItem = QGraphicsLineItem(0, self.coordY(prev), self.col_w-4,
                                               self.coordY(y), parent=self.item)
            lineItem.setX(x - self.col_w+2)
            lineItem.setPen(QPen(QColor(self.colors[i]),2))


