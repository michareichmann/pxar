#!/usr/bin/env python
# --------------------------------------------------------
#       Class for all the ROOT drawing stuff
# created on February 15th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import PyConfig
PyConfig.IgnoreCommandLineOptions = True
from ROOT import TGraphErrors, TGaxis, TLatex, TGraphAsymmErrors, TCanvas, TLegend, TArrow, TPad, TCutG, TLine, TPaveText, TPaveStats
from ROOT import gROOT, gStyle, kGreen, kOrange, kViolet, kYellow, kRed, kBlue, kMagenta, kAzure, kCyan, kTeal, TColor
from utils import *
from os.path import dirname, join
from numpy import ndarray, zeros, array
from uncertainties.core import Variable, ufloat, AffineScalarFunc
from screeninfo import get_monitors


class Draw:

    def __init__(self, verbose=True, titles=True, save_dir=''):

        # BASICS
        self.Verbose = verbose
        self.Res = load_resolution()
        self.ResultsDir = 'Results'
        self.SaveDir = save_dir
        self.FileTypes = ['root', 'png', 'pdf']

        # COLORS/SETTINGS
        self.count = 0
        self.colors = create_colorlist()
        self.FillColor = 821
        gStyle.SetLegendFont(42)
        self.Title = titles
        gStyle.SetOptTitle(titles)

        self.Objects = []

    # ----------------------------------------
    # region BASIC
    def set_save_directory(self, name):
        self.ResultsDir = name

    def make_bias_string(self, bias=None):
        if bias is None:
            return self.make_bias_string(self.bias) if hasattr(self, 'bias') else ''
        pol = 'm' if bias < 0 else 'p'
        return '_{pol}{bias:04d}'.format(pol=pol, bias=int(abs(bias)))

    def get_color(self):
        self.count %= 20
        color = self.colors[self.count]
        self.count += 1
        return color

    def get_colors(self, n):
        return array([self.get_color() for _ in xrange(n)], 'i')

    def reset_colors(self):
        self.count = 0
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region DRAWING
    def draw_histo(self, histo, save_name='', show=True, sub_dir=None, lm=.1, rm=.03, bm=.1, tm=None, draw_opt='', x=None, y=None, leg=None, logy=False, logx=False, logz=False,
                   canvas=None, grid=False, gridy=False, gridx=False, prnt=True, phi=None, theta=None):
        return self.save_histo(histo, save_name, show, sub_dir, lm, rm, bm, tm, draw_opt, x, y, leg, logy, logx, logz, canvas, grid, gridx, gridy, False, prnt, phi, theta)

    def draw_axis(self, x1, x2, y1, y2, title, limits=None, name='ax', col=1, width=1, off=.15, tit_size=.035, lab_size=0.035, tick_size=0.03, line=False, opt='+SU', l_off=.01, log=False):
        limits = ([y1, y2] if x1 == x2 else [x1, x2]) if limits is None else limits
        a = TGaxis(x1, y1, x2, y2, limits[0], limits[1], 510, opt + ('G' if log else ''))
        a.SetName(name)
        a.SetLineColor(col)
        a.SetLineWidth(width)
        a.SetLabelSize(lab_size if not line else 0)
        a.SetTitleSize(tit_size)
        a.SetTitleOffset(off)
        a.SetTitle(title)
        a.SetTitleColor(col)
        a.SetLabelColor(col)
        a.SetLabelFont(42)
        a.SetTitleFont(42)
        a.SetTickSize(tick_size if not line else 0)
        a.SetTickLength(tick_size if not line else 0)
        a.SetNdivisions(0) if line else do_nothing()
        a.SetLabelOffset(l_off)
        a.Draw()
        self.Objects.append(a)
        return a

    def draw_y_axis(self, x, ymin, ymax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return self.draw_axis(x, x, ymin, ymax, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    def draw_x_axis(self, y, xmin, xmax, tit, limits=None, name='ax', col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035, tick_size=0.03, l_off=.01, line=False, log=False):
        return self.draw_axis(xmin, xmax, y, y, tit, limits, name, col, w, off, tit_size, lab_size, tick_size, line, opt, l_off, log)

    def draw_line(self, x1, x2, y1, y2, color=1, width=1, style=1, name='li'):
        line = TCutG(name, 2, array([x1, x2], 'd'), array([y1, y2], 'd'))
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw('same')
        self.Objects.append(line)
        return line

    def draw_tline(self, x1, x2, y1, y2, color=1, width=1, style=1):
        line = TLine(x1, y1, x2, y2)
        line.SetLineColor(color)
        line.SetLineWidth(width)
        line.SetLineStyle(style)
        line.Draw()
        self.Objects.append(line)
        return line

    def draw_box(self, x1, y1, x2, y2, color=1, width=1, style=1, fillstyle=None, name='box', show=True):
        box = TCutG(name, 5, array([x1, x1, x2, x2, x1], 'd'), array([y1, y2, y2, y1, y1], 'd'))
        box.SetLineColor(color)
        box.SetFillColor(color)
        box.SetLineWidth(width)
        box.SetLineStyle(style)
        box.SetFillStyle(fillstyle) if fillstyle is not None else do_nothing()
        if show:
            box.Draw('same')
        self.Objects.append(box)
        return box

    def draw_vertical_line(self, x, ymin, ymax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(x, x, ymin, ymax, color, w, style, name) if not tline else self.draw_tline(x, x, ymin, ymax, color, w, style)

    def draw_horizontal_line(self, y, xmin, xmax, color=1, w=1, style=1, name='li', tline=False):
        return self.draw_line(xmin, xmax, y, y, color, w, style, name) if not tline else self.draw_tline(xmin, xmax, y, y, color, w, style)

    def draw_tlatex(self, x, y, text, name='text', align=20, color=1, size=.05, angle=None, ndc=None, font=None):
        latex = TLatex(x, y, text)
        self.format_text(latex, name, align, color, size, angle, ndc, font)
        latex.Draw()
        self.Objects.append(latex)
        return latex

    def draw_arrow(self, x1, x2, y1, y2, col=1, width=1, opt='<|', size=.005):
        ar = TArrow(x1, y1, x2, y2, size, opt)
        ar.SetLineWidth(width)
        ar.SetLineColor(col)
        ar.SetFillColor(col)
        ar.Draw()
        self.Objects.append(ar)

    def draw_tpad(self, name, tit='', pos=None, fill_col=0, gridx=None, gridy=None, margins=None, transparent=False, logy=None, logx=None, logz=None):
        margins = [.1, .1, .1, .1] if margins is None else margins
        pos = [0, 0, 1, 1] if pos is None else pos
        p = TPad(name, tit, *pos)
        p.SetFillColor(fill_col)
        p.SetMargin(*margins)
        do([p.SetLogx, p.SetLogy, p.SetLogz], [logx, logy, logz])
        do([p.SetGridx, p.SetGridy], [gridx, gridy])
        make_transparent(p) if transparent else do_nothing()
        p.Draw()
        p.cd()
        self.Objects.append(p)
        return p

    def draw_tpavetext(self, text, x1, x2, y1, y2, font=42, align=0, size=0, angle=0, margin=.05, color=1):
        p = TPaveText(x1, y1, x2, y2, 'ndc')
        p.SetFillColor(0)
        p.SetFillStyle(0)
        p.SetBorderSize(0)
        p.SetMargin(margin)
        t = p.AddText(text)
        self.format_text(t, 'pave', align, color, size, angle, ndc=True, font=font)
        p.Draw()
        self.Objects.append(p)
        return p

    def draw_preliminary(self, canvas=None, height=.06):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        return self.draw_tpavetext('#font[62]{RD42} Preliminary', c.GetLeftMargin(), .5, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=72, align=12, margin=0.04)

    def draw_irradiation(self, irr, canvas=None, height=.06, left=True):
        c = get_last_canvas() if canvas is None else canvas
        c.cd()
        x1, x2 = (c.GetLeftMargin(), .5) if left else (.5, 1 - c.GetRightMargin())
        return self.draw_tpavetext('Irradiation: {}'.format(irr), x1, x2, 1 - height - c.GetTopMargin(), 1 - c.GetTopMargin(), font=42, align=12, margin=0.04)

    def draw_stats(self, fit, y2=None, width=.3, prec='5.1f', names=None):
        names = fit.Names if names is None else names
        c = get_last_canvas()
        tm = .98 - .05 - c.GetTopMargin() if y2 is None else y2
        rm = .98 - c.GetRightMargin()
        p = TPaveStats(rm - width, tm - .06 * (fit.NPars + 1), rm, tm, 'ndc')
        p.SetBorderSize(1)
        p.SetFillColor(0)
        p.SetFillStyle(0)
        latex = p.AddText('Fit Result')
        latex.SetTextFont(42)
        ls = p.GetListOfLines()
        ls.Add(self.draw_tlatex(0, 0, '#chi^{{2}} / ndf  = {chi2:{p}} / {ndf}'.format(ndf=fit.Ndf(), chi2=fit.Chi2(), p=prec), size=0, align=0, font=42))
        for i in xrange(fit.NPars):
            ls.Add(self.draw_tlatex(0, 0, '{n}  = {v:{p}} #pm {e:{p}}'.format(n=names[i], v=fit.Parameter(i), e=fit.ParError(i), p=prec), size=0, align=0, font=42))
        p.Draw()
        self.Objects.append(p)
        return p

    def draw_frame(self, pad, xmin, xmax, ymin, ymax, tit, div=None, y_cent=None):
        pad.cd()
        fr = pad.DrawFrame(xmin, ymin, xmax, ymax)
        pad.Modified()
        fr.GetYaxis().SetTitle(tit)
        do(fr.GetYaxis().CenterTitle, y_cent)
        fr.GetYaxis().SetNdivisions(div) if div is not None else do_nothing()
        self.format_frame(fr)
        self.Objects.append(fr)
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region SAVING

    def save_plots(self, savename, sub_dir=None, canvas=None, prnt=True, save=True):
        """ Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
         such as in SSl terminals, it is recommended to pass the canvas to the method. """
        canvas = get_last_canvas() if canvas is None else canvas
        canvas.Modified()
        canvas.Update()
        if save:
            try:
                self.save_canvas(canvas, sub_dir=sub_dir, name=savename, print_names=prnt)
                self.Objects.append(canvas)
            except Exception as inst:
                warning('Error in save_canvas:\n{0}'.format(inst))

    def save_canvas(self, canvas, sub_dir=None, name=None, print_names=True):
        """ Saves the provided canvas into all the FileTypes. """
        sub_dir = self.SaveDir if sub_dir is None else sub_dir
        file_name = canvas.GetName() if name is None else name
        file_path = join(self.ResultsDir, sub_dir, '{typ}', file_name)
        set_root_output(False)
        for ext in self.FileTypes:
            ensure_dir(dirname(file_path.format(typ=ext)))
            canvas.SaveAs('{f}.{ext}'.format(f=file_path, ext=ext).format(typ=ext))
        if print_names:
            info('Saving plots: {nam}'.format(nam=file_name), prnt=self.Verbose)
        set_root_output(True)

    def save_histo(self, histo, save_name='test', show=True, sub_dir=None, lm=None, rm=None, bm=None, tm=None, draw_opt=None, x=None, y=None, leg=None, logy=False, logx=False, logz=False,
                   canvas=None, grid=False, gridx=False, gridy=False, save=True, prnt=True, phi=None, theta=None, sumw2=False):
        fac = 1 if self.Title else 1.16
        x = int(self.Res * fac) if x is None else int(x * self.Res)
        y = self.Res if y is None else int(y * self.Res)
        h = histo
        h.Sumw2(sumw2) if hasattr(h, 'Sumw2') and sumw2 is not None else do_nothing()
        set_root_output(show)
        c = TCanvas('c_{0}'.format(h.GetName()), h.GetTitle().split(';')[0], x, y) if canvas is None else canvas
        do(c.SetLeftMargin, lm)
        do(c.SetRightMargin, rm if rm is not None else None if round(c.GetRightMargin(), 1) != .1 else .03)
        do(c.SetBottomMargin, bm)
        do(c.SetTopMargin, tm if tm is not None else None if round(c.GetTopMargin(), 1) != .1 else .1 if self.Title else .03)
        c.SetLogx() if logx else do_nothing()
        c.SetLogy() if logy else do_nothing()
        c.SetLogz() if logz else do_nothing()
        c.SetGridx() if gridx or grid else do_nothing()
        c.SetGridy() if gridy or grid else do_nothing()
        c.SetPhi(phi) if phi is not None else do_nothing()
        c.SetTheta(theta) if theta is not None else do_nothing()
        h.Draw(draw_opt if draw_opt is not None else 'ap' if 'Graph' in h.ClassName() else '')
        if leg is not None:
            leg = [leg] if type(leg) is not list else leg
            for i in leg:
                i.Draw()
        self.save_plots(save_name, sub_dir=sub_dir, prnt=prnt, save=save)
        set_root_output(True)
        lst = [c, h, leg] if leg is not None else [c, h]
        self.Objects.append(lst)
        return c

    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region FORMATTING
    @staticmethod
    def format_text(t, name='text', align=20, color=1, size=.05, angle=None, ndc=None, font=None):
        t.SetName(name)
        t.SetTextAlign(align)
        t.SetTextColor(color)
        t.SetTextSize(size)
        do(t.SetTextAngle, angle)
        do(t.SetTextFont, font)
        do(t.SetNDC, ndc)
        return t

    @staticmethod
    def format_pie(pie, h=None, r=None, text_size=None, angle3d=None, angle_off=None, label_format=None):
        do([pie.SetHeight, pie.SetRadius], [h, r])
        do(pie.SetTextSize, text_size)
        do(pie.SetAngle3D, angle3d)
        do(pie.SetLabelFormat, label_format)
        do(pie.SetAngularOffset, angle_off)

    @staticmethod
    def format_frame(frame):
        fr = frame
        fr.GetYaxis().SetTitleSize(.06)
        fr.GetYaxis().SetTitleOffset(.6)
        fr.GetYaxis().SetLabelSize(.06)
        fr.SetTitleSize(.05)
        fr.GetXaxis().SetTickLength(0)
        fr.GetXaxis().SetLabelOffset(99)
        fr.SetLineColor(0)
        fr.GetXaxis().SetTimeDisplay(1)
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region MAKING
    def make_tgrapherrors(self, name, title, color=1, marker=20, marker_size=1, width=1, asym_err=False, style=1, x=None, y=None, ex=None, ey=None):
        x = list(x) if type(x) == ndarray else x
        if x is None:
            gr = TGraphErrors() if not asym_err else TGraphAsymmErrors()
        else:
            gr = TGraphErrors(*make_graph_args(x, y, ex, ey)) if not asym_err else TGraphAsymmErrors(*make_graph_args(x, y, ex, ey))
        gr.SetTitle(title)
        gr.SetName(name)
        gr.SetMarkerStyle(marker)
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)
        gr.SetMarkerSize(marker_size)
        gr.SetLineWidth(width)
        gr.SetLineStyle(style)
        self.Objects.append(gr)
        return gr

    def make_legend(self, x1=.65, y2=.88, nentries=2, scale=1, name='l', y1=None, clean=False, margin=.25, x2=None):
        x2 = .95 if x2 is None else x2
        y1 = y2 - nentries * .05 * scale if y1 is None else y1
        legend = TLegend(x1, y1, x2, y2)
        legend.SetName(name)
        legend.SetTextFont(42)
        legend.SetTextSize(0.03 * scale)
        legend.SetMargin(margin)
        if clean:
            legend.SetLineWidth(2)
            legend.SetBorderSize(0)
            legend.SetFillColor(0)
            legend.SetFillStyle(0)
            legend.SetTextAlign(12)
        self.Objects.append(legend)
        return legend

    def make_canvas(self, name='c', title='c', x=1., y=1., show=True, logx=None, logy=None, logz=None, gridx=None, gridy=None, transp=None):
        set_root_output(show)
        c = TCanvas(name, title, int(x * self.Res), int(y * self.Res))
        do([c.SetLogx, c.SetLogy, c.SetLogz], [logx, logy, logz])
        do([c.SetGridx, c.SetGridy], [gridx, gridy])
        do(make_transparent, c, transp)
        self.Objects.append(c)
        return c

    def make_graph_from_profile(self, p):
        x_range = [i for i in xrange(p.GetNbinsX()) if p.GetBinContent(i)]
        x = [make_ufloat([p.GetBinCenter(i), p.GetBinWidth(i) / 2]) for i in x_range]
        y = [make_ufloat([p.GetBinContent(i), p.GetBinError(i)]) for i in x_range]
        return self.make_tgrapherrors('g{n}'.format(n=p.GetName()[1:]), p.GetTitle(), x=x, y=y)
    # endregion
    # ----------------------------------------

    def format_statbox(self, x=.95, y=None, w=.2, n_entries=3, only_fit=False, fit=False, entries=False, form=None, m=False, rms=False, all_stat=False):
        gStyle.SetOptFit(int(only_fit or fit))
        opt_stat = '100000{}{}{}0'.format(*[1 if val else 0 for val in [rms, m, entries]] if not all_stat else [1, 1, 1])
        if only_fit:
            opt_stat = '0011'
        if fit:
            opt_stat = '1111'
        y = (.88 if self.Title else .95) if y is None else y
        gStyle.SetOptStat(int(opt_stat))
        gStyle.SetFitFormat(form) if form is not None else do_nothing()
        gStyle.SetStatX(x)
        gStyle.SetStatY(y)
        gStyle.SetStatW(w)
        gStyle.SetStatH(.04 * n_entries)


def create_colorlist():
    col_names = [kGreen, kOrange, kViolet, kYellow, kRed, kBlue, kMagenta, kAzure, kCyan, kTeal]
    colors = []
    for color in col_names:
        colors.append(color + (1 if color != 632 else -7))
    for color in col_names:
        colors.append(color + (3 if color != 800 else 9))
    return colors


def make_graph_args(x, y, ex=None, ey=None):
    if len(list(x)) != len(list(y)):
        warning('Arrays have different size!')
        return []
    lx = len(x)
    x = x if type(x[0]) in [Variable, AffineScalarFunc] else [make_ufloat(tup) for tup in zip(x, zeros(lx) if ex is None else ex)]
    y = y if type(y[0]) in [Variable, AffineScalarFunc] else [make_ufloat(tup) for tup in zip(y, zeros(lx) if ey is None else ey)]
    return [lx, array([v.n for v in x], 'd'), array([v.n for v in y], 'd'), array([v.s for v in x], 'd'), array([v.s for v in y], 'd')]


def make_transparent(pad):
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.SetFrameFillStyle(4000)


def get_last_canvas():
    try:
        return gROOT.GetListOfCanvases()[-1]
    except IndexError:
        warning('There is no canvas is in the list...')


def set_root_output(status=True):
    gROOT.SetBatch(not status)
    gROOT.ProcessLine('gErrorIgnoreLevel = {e};'.format(e='0' if status else 'kError'))


def set_time_axis(histo, form='%H:%M', off=0):
    histo.GetXaxis().SetTimeFormat(form)
    histo.GetXaxis().SetTimeOffset(-off - 3600 if off else 0)
    histo.GetXaxis().SetTimeDisplay(1)


def set_palette(pal=1, custom=False):
    if custom:
        stops = array([0., .5, 1], 'd')
        green = array([0. / 255., 200. / 255., 80. / 255.], 'd')
        blue = array([0. / 255., 0. / 255., 0. / 255.], 'd')
        red = array([180. / 255., 200. / 255., 0. / 255.], 'd')
        gStyle.SetNumberContours(20)
        bla = TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, 255)
        color_table = array([bla + ij for ij in xrange(255)], 'i')
        gStyle.SetPalette(len(color_table), color_table)
    else:
        gStyle.SetPalette(pal)


def load_resolution():
    try:
        m = get_monitors()
        return round_down_to(m[0].height, 500)
    except Exception as err:
        warning(err)
        return 1000


def make_ufloat(tup):
    if type(tup) in [Variable, AffineScalarFunc]:
        return tup
    if type(tup) in [tuple, list, ndarray]:
        return ufloat(*tup) if type(tup[0]) not in [Variable, AffineScalarFunc] else ufloat(tup[0].n, tup[1].n)
    return ufloat(tup, 0)


def set_z_range(zmin, zmax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetZaxis().SetRangeUser(zmin, zmax)


def format_histo(histo, name=None, title=None, x_tit=None, y_tit=None, z_tit=None, marker=20, color=None, line_color=None, markersize=None, x_off=None, y_off=None, z_off=None, lw=1,
                 fill_color=None, fill_style=None, stats=True, tit_size=None, lab_size=None, l_off_y=None, l_off_x=None, draw_first=False, x_range=None, y_range=None, z_range=None, sumw2=None,
                 do_marker=True, style=None, ndivx=None, ndivy=None, ncont=None, tick_size=None, t_ax_off=None, center_y=False, center_x=False, yax_col=None, normalise=None, pal=None, rebin=None):
    h = histo
    if draw_first:
        set_root_output(False)
        h.Draw('nostack' if h.ClassName() == 'THStack' else 'a')
        set_root_output(True)
    do(h.SetTitle, title)
    do(h.SetName, name)
    do(set_palette, pal)
    if normalise is not None:
        y_tit = y_tit.replace('Number', 'Percentage') if y_tit is not None else y_tit
        h.Sumw2(True)
        normalise_histo(h)
    try:
        h.SetStats(stats)
    except AttributeError or ReferenceError:
        pass
    do(h.Rebin, rebin) if hasattr(h, 'Rebin') else do_nothing()
    # markers
    try:
        if do_marker:
            do(h.SetMarkerStyle, marker)
            do(h.SetMarkerColor, color)
            do(h.SetMarkerSize, markersize)
    except AttributeError or ReferenceError:
        pass
    # lines/fill
    try:
        h.SetLineColor(line_color) if line_color is not None else h.SetLineColor(color) if color is not None else do_nothing()
        h.SetLineWidth(lw)
        h.SetFillColor(fill_color) if fill_color is not None else do_nothing()
        h.SetFillStyle(fill_style) if fill_style is not None else do_nothing()
        h.SetFillStyle(style) if style is not None else do_nothing()
        h.SetContour(ncont) if ncont is not None else do_nothing()
    except AttributeError or ReferenceError:
        pass
    # axes
    try:
        x_args = [x_tit, x_off, tit_size, center_x, lab_size, l_off_x, x_range, ndivx, tick_size, ]
        y_args = [y_tit, y_off, tit_size, center_y, lab_size, l_off_y, y_range, ndivy, tick_size, yax_col]
        z_args = [z_tit, z_off, tit_size, False, lab_size, None, z_range, None, tick_size]
        for i, name in enumerate(['X', 'Y', 'Z']):
            format_axis(getattr(h, 'Get{}axis'.format(name))(), is_graph(h), *[x_args, y_args, z_args][i])
    except AttributeError or ReferenceError:
        pass
    set_time_axis(h, off=t_ax_off) if t_ax_off is not None else do_nothing()
    do(h.Sumw2, sumw2) if hasattr(h, 'Sumw2') else do_nothing()


def format_axis(axis, graph, title, tit_offset, tit_size, centre_title, lab_size, label_offset, limits, ndiv, tick_size, color=None):
    do(axis.SetTitle, title)
    do(axis.SetTitleOffset, tit_offset)
    do(axis.SetTitleSize, tit_size)
    axis.CenterTitle(centre_title)
    do(axis.SetLabelSize, lab_size)
    do(axis.SetLabelOffset, label_offset)
    if limits is not None:
        axis.SetLimits(*limits) if graph and 'xaxis' in axis.GetName() else axis.SetRangeUser(*limits)
    do(axis.SetNdivisions, ndiv)
    do(axis.SetTickSize, tick_size)
    do(axis.SetTitleColor, color)
    do(axis.SetLabelColor, color)
    do(axis.SetAxisColor, color)


def normalise_histo(histo, x_range=None, from_min=False):
    h = histo
    x_axis = h.GetXaxis()
    x_axis.SetRangeUser(*x_range) if x_range is not None else do_nothing()
    min_bin = h.GetMinimumBin() if from_min else 0
    integral = h.Integral(min_bin, h.GetNbinsX() - 1)
    return scale_histo(h, integral)


def scale_histo(histo, value=None, to_max=False, x_range=None):
    h = histo
    maximum = h.GetBinContent(h.GetMaximumBin())
    if x_range is not None:
        h.GetXaxis().SetRangeUser(*x_range) if x_range is not None else do_nothing()
        maximum = h.GetBinContent(h.GetMaximumBin())
        h.GetXaxis().UnZoom()
    value = maximum if to_max else value
    if value:
        h.Scale(1. / value)
    return h


def is_graph(h):
    return 'Graph' in h.ClassName()
