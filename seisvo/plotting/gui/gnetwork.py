import os
import datetime as dt
from seisvo import DB_PATH
from seisvo.core.obspyext import Stream2

class NetworkCanvas(FigureCanvas):
    def __init__(self, network, sta_list, starttime, delta, component="Z", sde_file=None, parent=None, **kwargs):

        # load objects
        self.network = network
        self.sta_list = sta_list
        self.starttime = starttime
        self.component = component
        self.delta = delta
        self.parent = parent # mainwindow of the GUI

        # self.endtime = starttime + dt.timedelta(minutes=delta)

        # kwargs        
        self.olap = kwargs.get("delta_olap", 0.2)
        self.samples = kwargs.get("samples", 50)
        
        self.plot_kwargs = dict(
            interval = kwargs.get("interval", 60),
            fq_band = kwargs.get("fq_band", [0.5, 10]),
            color_name = kwargs.get("color_name", 'zesty'),
            remove_response = kwargs.get("remove_response", False),
            sample_rate = kwargs.get("sample_rate", 40),
        )

        # load specgram
        specgram = kwargs.get("specgram", None)
        if specgram:
            if specgram not in [st.stats.id for st in station_list]:
                print(" warn: no specgram loaded")
                specgram = None
        self.plot_kwargs["specgram"] = specgram

        self.psd_frame = None

        self.click_info_dict = {
            'left':None,
            'right':None,
            'diff':None,
            'freq':None,
            'trace':None
        }

        if not sde_file:
            sde_file = os.path.join(DB_PATH, network.stats.code) + '.db'
        try:
            self.sde = SDE(sde_file)
        except:
            raise ValueError(f" error loading SDE {sde_file}. Revise file/location")
        
        # load canvas
        self.fig = Figure(figsize=(20,9))
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, 
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self.mpl_connect('button_release_event', self.on_click)
        self.mpl_connect('key_press_event', self.on_key)
        self.setPlot()
    

    def setPlot(self):
        self.fig.clf()
        self.parent.plainTextEdit.setText('')

        with pyqtgraph.BusyCursor():
            ans = get_fig(
                self.station, 
                self.starttime,
                self.channel,
                self.starttime + dt.timedelta(minutes=delta),
                return_axes=True,
                fig=self.fig,
                **self.plot_kwargs
                )
        
        self.chan_axes, self.specgram_axes, self.time, self.trace_list = ans

        if self.specgram_axes:
            self.nav = Navigation(self.chan_axes, imshow_axes=self.specgram_axes[0], parent=self)
        else:
            self.nav = Navigation(self.chan_axes, parent=self)

        self.__eventlist = {}
        self.__show_events()

        self.showCatalog = False
        with pyqtgraph.BusyCursor():
            self.draw()


    def showTickInfo(self):
        text = f" {self.info_dict['trace']}:\n"
        text += f" L: {self.info_dict['left']}\n"
        text += f" R: {self.info_dict['right']}\n"

        if self.info_dict['diff']:
            text += f" R-L: {self.info_dict['diff']:.2f} sec\n"
        
        if self.info_dict['freq']:
            text += f" Freq.: {self.info_dict['freq']:.2f}\n"
        self.parent.plainTextEdit.setText(text)

    


def get_text(remove_response, x):
    if remove_response:
        txt = f'Max.: {x*10**6:.2f}'
        return txt + r' $\mu m$'
    else:
        txt = f'Max.: {x:.2f}'
        return txt + ' cnts'


def get_fig(station_list, starttime, component, endtime, return_axes=False, **kwargs):
    """[summary]

    Args:
        station: List of station objects
        starttime: datetime
        channel: the channel/s to plot
        endtime: datetime
    """

    # define some variables
    fq_band = kwargs.get('fq_band')
    remove_response = kwargs.get('remove_response')
    sample_rate = kwargs.get('sample_rate')
    color_name = kwargs.get('color_name')
    specgram = kwargs.get('specgram')

    colors = get_colors(color_name)

    # define the stream
    stream = Stream2()
    for n, sta in enumerate(station_list):
        if component in ["Z", "N", "E"]:
            channel = sta.get_chan(component)
        else:
            raise ValueError("multiple channel is under development: use gstation instead")
        
        try:
            stream += sta.get_stream(starttime, endtime, sample_rate=sample_rate, chan=channel, remove_response=remove_response)
        except:
            print(f" warn: error loading station {sta.stats.id}")

    # define the figure frame
    fig = kwargs.get('fig', None)
    grid = {'hspace':0.1, 'left':0.08, 'right':0.92, 'wspace':0.05, 'top':0.95, 'bottom':0.05}
    
    nrows = len(stream)
    ncols = 1

    if specgram in ['.'.join(s.id.split('.')[1:3]) for s in stream]:
        nrows += 1
        ncols += 1
        grid['width_ratios'] = [1, 0.01]

    if not fig:
        dpi = kwargs.get('dpi', 100)
        figsize = kwargs.get('figsize', (20,9))
        fig, axes = plt.subplots(nrows, ncols, gridspec_kw=grid, figsize=figsize, dpi=dpi)
    else:
        axes = fig.subplots(nrows, ncols, gridspec_kw=grid)
    
    if isinstance(axes, np.ndarray):
        axes = axes.reshape(nrows, ncols)
    else:
        axes = np.array([axes]).reshape(nrows, ncols)
    
    # if specgram, plot first
    i = 0
    spec_return = ()
    specgram_trace = stream.get_component(component, station=specgram.split('.')[1], loc=specgram.split('.')[2])
    if specgram_trace:
        if remove_response:
            label = r"[dB cnts$^2$/Hz]"
        else:
            label = r"[dB m$^2$s$^{-2}$/Hz]"
        spec_ax = axes[i,0]
        specgram_trace.specgram(
            axes=spec_ax, 
            fq_band=fq_band, 
            per_lap=0.75,
            xlabel=' ',
            axis_bar=spec_bar_ax, 
            axis_bar_label=label
            )
        spec_ax.set_ylabel(spectram, fontsize=12, color=colors[specgram])
        spec_ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        spec_ax.xaxis.set_major_formatter(mtick.NullFormatter())
        spec_return = (spec_ax, spec_bar_ax)
        i += 1
    else:
        specgram_trace = stream[0]

    # plot channels, being the first channel the same used for specgram
    chan_axes = []
    time = None # read time once
    trace_list = []
    while i < len(stream):
        for trace in stream:
            sta_id = '.'.join(trace.stats.network, trace.stats.station, trace.stats.location)
            txt = '.'.join(trace.stats.station, trace.stats.location, trace.stats.channel[-1])
            trace_list += [txt]
            color_tr = colors[sta_id]
            if trace.id == specgram_trace.id:
                ax = axes[i,0]
                i += 1
            else:
                ax = axes[i,0]
                i += 1

            chan_axes += [ax]
            
            if ncols == 2:
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)
                ax.set_frame_on(False)
            
            if time == None:
                time = trace.get_time()

            data = trace.get_data(detrend=True, fq_band=fq_band)
            max_ampl = np.nanmax(np.abs(data))
            norm_data = data/max_ampl
            
            ax.plot(time, norm_data, color=color_tr, lw=1.0)
            ax.set_ylabel(txt, fontsize=12, color=color_tr)
            ax.set_ylim(-1, 1)

            ax.annotate(
                get_text(remove_response, max_ampl),
                xy=(-0.01,0.75),
                xycoords='axes fraction',
                color=color_tr,
                bbox=dict(boxstyle="round", fc="w", alpha=1)
            )

            if i == len(stream):
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            else:
                ax.xaxis.set_major_formatter(mtick.NullFormatter())

            ax.yaxis.set_major_locator(mtick.NullLocator())
            ax.grid(axis='x', which='major', color='k', ls='--', alpha=0.4)
            ax.grid(axis='x', which='minor', color='k', ls='--', alpha=0.2)

    if spec_return:
        xaxis = [spec_ax] + chan_axes
    else:
        xaxis = chan_axes
    
    [ax.set_xlim(time[0], time[-1]) for ax in xaxis]
    [ax.xaxis.set_major_locator(mtick.MaxNLocator(nbins=6)) for ax in xaxis]
    [ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(4)) for ax in xaxis]

    if (time[-1] - time[0]).total_seconds() > 86400:
        xaxis[0].set_title('%s--%s', time[0].strftime('%d %b %Y'), time[-1].strftime('%d %b %Y'))
    else:
        xaxis[0].set_title(time[0].strftime('%d %b %Y'))
    
    fig.align_ylabels()

    if return_axes:
        return chan_axes, spec_return, time, trace_list

    else:
        return fig, axes
        

        

