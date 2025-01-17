# distutils: language = c++
from libcpp cimport bool
from libc.stdint cimport uint8_t, int8_t, uint16_t, int16_t, int32_t, uint32_t
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.map cimport map
import numpy

cimport PyPxarCore

FLAG_FORCE_SERIAL   = int(_flag_force_serial)
FLAG_CALS           = int(_flag_cals)
FLAG_XTALK          = int(_flag_xtalk)
FLAG_RISING_EDGE    = int(_flag_rising_edge)
FLAG_DISABLE_DACCAL = int(_flag_disable_daccal)
FLAG_NOSORT         = int(_flag_nosort)
FLAG_CHECK_ORDER    = int(_flag_check_order)
FLAG_FORCE_UNMASKED = int(_flag_force_unmasked)
FLAG_DUMP_FLAWED_EVENTS = int(_flag_dump_flawed_events)
FLAG_DISABLE_READBACK_COLLECTION = int(_flag_disable_readback_collection)
FLAG_DISABLE_EVENTID_CHECK = int(_flag_disable_eventid_check)
FLAG_ENABLE_XORSUM_LOGGING = int(_flag_enable_xorsum_logging)

cdef class Pixel:
    cdef pixel *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, address = None, data = None): # default to None to mimick overloading of constructor
        if address is not None and data is not None:
            self.thisptr = new pixel(address, data)
        else:
            self.thisptr = new pixel()
    def __dealloc__(self):
        del self.thisptr
    def __str__(self):
        s = "ROC " + str(self.roc)
        s += " [" + str(self.column) + "," + str(self.row) + "," + str(self.value) + "] "
        return s
    cdef fill(self, pixel p):
        self.thisptr.setRoc(p.roc())
        self.thisptr.setColumn(p.column())
        self.thisptr.setRow(p.row())
        self.thisptr.setValue(p.value())
        self.thisptr.setBufferCorruption(p.bufferCorruption())
        self.thisptr.setInvalidAddress(p.invalidAddress())
        self.thisptr.setInvalidPulseHeight(p.invalidPulseHeight())
    cdef c_clone(self, pixel* p):
        del self.thisptr
        self.thisptr = p
    property roc:
        def __get__(self): return self.thisptr.roc()
        def __set__(self, roc): self.thisptr.setRoc(roc)
    property column:
        def __get__(self): return self.thisptr.column()
        def __set__(self, column): self.thisptr.setColumn(column)
    property row:
        def __get__(self): return self.thisptr.row()
        def __set__(self, row): self.thisptr.setRow(row)
    property value:
        def __get__(self): return self.thisptr.value()
        def __set__(self, value): self.thisptr.setValue(value)
    property buffer_corruption:
        def __get__(self): return self.thisptr.bufferCorruption()
        def __set__(self, value): self.thisptr.setBufferCorruption(value)
    property invalid_address:
        def __get__(self): return self.thisptr.invalidAddress()
        def __set__(self, value): self.thisptr.setInvalidAddress(value)
    property invalid_pulse_height:
        def __get__(self): return self.thisptr.invalidPulseHeight()
        def __set__(self, value): self.thisptr.setInvalidPulseHeight(value)

cdef class PixelConfig:
    cdef pixelConfig *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, column = None, row = None, trim = None): # default to None to mimick overloading of constructor
        if column is not None and row is not None and trim is not None:
            self.thisptr = new pixelConfig(column, row, trim)
        else:
            self.thisptr = new pixelConfig()
    def __dealloc__(self):
        del self.thisptr
    def __str__(self):
        s = "ROC " + str(self.roc) + " [" + str(self.column) + "," + str(self.row) + "] Trim: " + str(self.trim) + " Mask: "
        s += "True" if self.mask else "False"
        return s
    def __richcmp__(self, other not None, int op):
        if op == 2: # ==
            return (self.roc == other.roc and self.column == other.column and self.row == other.row)
        elif op == 3: # !=
            return (self.roc != other.roc or self.column != other.column or self.row != other.row)
    cdef fill(self, pixelConfig p):
        self.thisptr.setRoc(p.roc())
        self.thisptr.setColumn(p.column())
        self.thisptr.setRow(p.row())
        self.thisptr.setTrim(p.trim())
        self.thisptr.setEnable(p.enable())
        self.thisptr.setMask(p.mask())
    cdef c_clone(self, pixelConfig* p):
        del self.thisptr
        thisptr = p
    property roc:
        def __get__(self): return self.thisptr.roc()
        def __set__(self, roc): self.thisptr.setRoc(roc)
    property column:
        def __get__(self): return self.thisptr.column()
        def __set__(self, column): self.thisptr.setColumn(column)
    property row:
        def __get__(self): return self.thisptr.row()
        def __set__(self, row): self.thisptr.setRow(row)
    property trim:
        def __get__(self): return self.thisptr.trim()
        def __set__(self, trim): self.thisptr.setTrim(trim)
    property enable:
        def __get__(self): return self.thisptr.enable()
        def __set__(self, enable): self.thisptr.setEnable(enable)
    property mask:
        def __get__(self): return self.thisptr.mask()
        def __set__(self, mask): self.thisptr.setMask(mask)

    def __repr__(self):
        return self.__str__()


cdef class RocConfig:
    cdef rocConfig *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
            self.thisptr = new rocConfig()
    def __dealloc__(self):
        del self.thisptr
    property column:
        def __get__(self):
            r = list()
            for p in self.thisptr.pixels:
                P = PixelConfig()
                P.c_clone(&p)
                r.append(P)
            return r
        def __set__(self, value):
            cdef vector[pixelConfig] v
            cdef PixelConfig pc
            for pc in value:
                v.push_back( <pixelConfig> pc.thisptr[0])
            self.thisptr.pixels = v
    property dacs:
        def __get__(self): return self.thisptr.dacs
        def __set__(self, value):
            cdef map[uint8_t, uint8_t] m
            for key in value.iterkeys():
                m[key] = value[key]
            self.thisptr.dacs = m
    property type:
        def __get__(self): return self.thisptr.type
        def __set__(self, value): self.thisptr.type = value
    property enable:
        def __get__(self): return self.thisptr.enable()
        def __set__(self, enable): self.thisptr.setEnable(enable)

cdef class Statistics:
    cdef statistics thisobj      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisobj = statistics()
    def __dealloc__(self):
        pass
    def __str__(self):
        return "Decoding statistics..."
    cdef c_clone(self, statistics s):
        self.thisobj = s
    def dump(self):
        self.thisobj.dump()
    property errors:
        def __get__(self): return self.thisobj.errors()
    property info_words_read:
        def __get__(self): return self.thisobj.info_words_read()
    property errors_event:
        def __get__(self): return self.thisobj.errors_event()
    property errors_tbm:
        def __get__(self): return self.thisobj.errors_tbm()
    property errors_roc:
        def __get__(self): return self.thisobj.errors_roc()
    property errors_pixel:
        def __get__(self): return self.thisobj.errors_pixel()
    property valid_pixels:
        def __get__(self): return self.thisobj.info_pixels_valid()
    property total_events:
        def __get__(self): return self.thisobj.info_events_total()
    property valid_events:
        def __get__(self): return self.thisobj.info_events_valid()
    property empty_events:
        def __get__(self): return self.thisobj.info_events_empty()
    property errors_event_start:
        def __get__(self): return self.thisobj.errors_event_start()
    property errors_event_stop:
        def __get__(self): return self.thisobj.errors_event_stop()
    property errors_event_overflow:
        def __get__(self): return self.thisobj.errors_event_overflow()
    property errors_event_invalid_words:
        def __get__(self): return self.thisobj.errors_event_invalid_words()
    property errors_event_invalid_xor:
        def __get__(self): return self.thisobj.errors_event_invalid_xor()
    property errors_event_frame:
        def __get__(self): return self.thisobj.errors_event_frame()
    property errors_event_idledata:
        def __get__(self): return self.thisobj.errors_event_idledata()
    property errors_event_nodata:
        def __get__(self): return self.thisobj.errors_event_nodata()
    property errors_event_pkam:
        def __get__(self): return self.thisobj.errors_event_pkam()
    property errors_tbm_header:
        def __get__(self): return self.thisobj.errors_tbm_header()
    property errors_tbm_eventid_mismatch:
        def __get__(self): return self.thisobj.errors_tbm_eventid_mismatch()
    property errors_tbm_trailer:
        def __get__(self): return self.thisobj.errors_tbm_trailer()
    property errors_roc_missing:
        def __get__(self): return self.thisobj.errors_roc_missing()
    property errors_roc_readback:
        def __get__(self): return self.thisobj.errors_roc_readback()
    property errors_pixel_incomplete:
        def __get__(self): return self.thisobj.errors_pixel_incomplete()
    property errors_pixel_address:
        def __get__(self): return self.thisobj.errors_pixel_address()
    property errors_pixel_pulseheight:
        def __get__(self): return self.thisobj.errors_pixel_pulseheight()
    property errors_pixel_buffer_corrupt:
        def __get__(self): return self.thisobj.errors_pixel_buffer_corrupt()

cdef class PxEvent:
    cdef Event *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new Event()
    def __dealloc__(self):
        del self.thisptr
    def __str__(self):
        s = "====== "
        for i in self.header: s += hex(i) + " "
        s += " ====== "
        for px in self.pixels: s += str(px)
        s += " ====== "
        for i in self.trailer: s += hex(i) + " "
        s += " ======\n"
        return str(s)
    def __repr__(self):
        return '[{}]'.format(', '.join(px for px in self.pixels))
    cdef clone(self, Event ev):
        self.thisptr = new Event(ev)
    def printHeader(self):
        self.thisptr.printHeader()
    def printTrailer(self):
        self.thisptr.printTrailer()

    property pixels:
        def __get__(self):
            r = list()
            for p in self.thisptr.pixels:
                P = Pixel()
                P.fill(p)
                r.append(P)
            return r
        def __set__(self, value):
            cdef vector[pixel] v
            cdef Pixel px
            for px in value:
                v.push_back( <pixel> px.thisptr[0])
            self.thisptr.pixels = v
    property header:
        def __get__(self): return self.thisptr.getHeaders()
        def __set__(self, value): self.thisptr.addHeader(value)
    property trailer:
        def __get__(self): return self.thisptr.getTrailers()
        def __set__(self, trailer): self.thisptr.addTrailer(trailer)
    property haveNoTokenPass:
        def __get__(self): return self.thisptr.haveNoTokenPass()
    property haveTokenPass:
        def __get__(self): return self.thisptr.haveTokenPass()
    property haveResetTBM:
        def __get__(self): return self.thisptr.haveResetTBM()
    property haveResetROC:
        def __get__(self): return self.thisptr.haveResetROC()
    property haveSyncError:
        def __get__(self): return self.thisptr.haveSyncError()
    property haveSyncTrigger:
        def __get__(self): return self.thisptr.haveSyncTrigger()
    property haveClearTriggerCount:
        def __get__(self): return self.thisptr.haveClearTriggerCount()
    property haveCalTrigger:
        def __get__(self): return self.thisptr.haveCalTrigger()
    property stacksFull:
        def __get__(self): return self.thisptr.stacksFull()
    property haveAutoReset:
        def __get__(self): return self.thisptr.haveAutoReset()
    property havePkamReset:
        def __get__(self): return self.thisptr.havePkamReset()
    property triggerCounts:
        def __get__(self): return self.thisptr.triggerCounts()
    property triggerPhases:
        def __get__(self): return self.thisptr.triggerPhases()
    property dataIDs:
        def __get__(self): return self.thisptr.dataIDs()
    property dataValues:
        def __get__(self): return self.thisptr.dataValues()
    property stackCounts:
        def __get__(self): return self.thisptr.stackCounts()
    property incomplete_data:
        def __get__(self): return self.thisptr.incomplete_data
    property missing_roc_headers:
        def __get__(self): return self.thisptr.missing_roc_headers
    property roc_readback:
        def __get__(self): return self.thisptr.roc_readback
    property no_data:
        def __get__(self): return self.thisptr.no_data
    property eventid_mismatch:
        def __get__(self): return self.thisptr.eventid_mismatch

cdef class PyPxarCore:
    cdef pxarCore *thisptr # hold the C++ instance
    def __cinit__(self, usbId = "*", logLevel = "INFO"):
        self.thisptr = new pxarCore(usbId, logLevel)
    def __dealloc__(self):
        del self.thisptr
    def initTestboard(self,sig_delays, power_settings, pg_setup):
        """ Initializer method for the testboard
        Parameters are dictionaries in the form {"name":value}:
        sig_delays = signal delays
        power_settings = voltage/current limit settings
        pg_setup = initial pattern generator setup
        """
        cdef vector[pair[string, uint8_t]] sd
        cdef vector[pair[string, double]] ps
        cdef vector[pair[string, uint8_t ]] pgs
        # type conversions for fixed-width integers need to
        # be handled very explicitly: creating pairs to push into vects
        for key, value in sig_delays.items():
            sd.push_back(pair[string,uint8_t](key,int(value)))
        for key, value in power_settings.items():
            ps.push_back((key,float(value)))
        for key, value in pg_setup:
            pgs.push_back(pair[string, uint8_t](key, int(value)))
        return self.thisptr.initTestboard(sd, ps, pgs)
    def setTestboardPower(self, power_settings):
        """ Initializer method for the testboard
        Parameters are dictionaries in the form {"name":value}:
        power_settings = voltage/current limit settings
        """
        cdef vector[pair[string, double]] ps
        # type conversions for fixed-width integers need to
        # be handled very explicitly: creating pairs to push into vects
        for key, value in power_settings.items():
            ps.push_back((key,float(value)))
        self.thisptr.setTestboardPower(ps)
    def setBlackOffsets(self, values):
        cdef vector[float] offsets
        for offset in values:
            offsets.push_back(offset)
        self.thisptr.setBlackOffsets(offsets)
    def setDecodingL1Offsets(self, values):
        cdef vector[float] offsets
        for offset in values:
            offsets.push_back(offset)
        self.thisptr.setDecodingL1Offsets(offsets)
    def setDecodingAlphas(self, values):
        cdef vector[float] alphas
        for alpha in values:
            alphas.push_back(alpha)
        self.thisptr.setDecodingAlphas(alphas)
    def getTestboardDelays(self):
        r = self.thisptr.getTestboardDelays()
        return {tup.first: tup.second for tup in r}
    def setTestboardDelays(self, sig_delays):
        """ Initializer method for the testboard
        Parameters are dictionaries in the form {"name":value}:
        sig_delays = signal delays
        """
        cdef vector[pair[string, uint8_t]] sd
        # type conversions for fixed-width integers need to
        # be handled very explicitly: creating pairs to push into vects
        for key, value in sig_delays.items():
            sd.push_back(pair[string,uint8_t](key,value))
        self.thisptr.setTestboardDelays(sd)
    def setPatternGenerator(self, pg_setup):
        """ Initializer method for the testboard
        Parameters are dictionaries in the form {"name":value}:
        pg_setup = initial pattern generator setup
        """
        cdef vector[pair[string, uint8_t ]] pgs
        # type conversions for fixed-width integers need to
        # be handled very explicitly: creating pairs to push into vects
        for item in enumerate(pg_setup):
            pgs.push_back(pair[string, uint8_t ](item[1][0],item[1][1]))
        self.thisptr.setPatternGenerator(pgs)
    def initDUT(self, hubids, tbmtype, tbmDACs, roctype, rocDACs, rocPixels, rocI2C = None):
        """ Initializer method for the DUT (attached devices)
        Parameters:
	hubId (int vector)
        tbmtype (string)
        tbmDACs (list of dictionaries (string,int), one for each TBM)
        roctype (string)
        rocDACs (list of dictionaries (string,int), one for each ROC)
        rocPixels (list of list of pixelConfigs, one list for each ROC)
        rocI2C (list of I2C addresses of the ROCs)
        """
        cdef vector[uint8_t] hubs
        cdef vector[vector[pair[string,uint8_t]]] td
        cdef vector[vector[pair[string,uint8_t]]] rd
        cdef vector[vector[pixelConfig]] rpcs
        cdef PixelConfig pc
        cdef vector[uint8_t] i2c

        if isinstance(hubids,list):
            for i in hubids:
                hubs.push_back(i)
        else:
            hubs.push_back(hubids)

        for idx, tbmDAC in enumerate(tbmDACs):
            td.push_back(vector[pair[string,uint8_t]]())
            for key, value in tbmDAC.items():
                td.at(idx).push_back(pair[string,uint8_t](key,int(value,0) if isinstance(value,str) else int(value)))
        for idx, rocDAC in enumerate(rocDACs):
            rd.push_back(vector[pair[string,uint8_t]]())
            for key, value in rocDAC.items():
                rd.at(idx).push_back(pair[string,uint8_t](key,int(value,0) if isinstance(value,str) else int(value)))
        for idx, rocPixel in enumerate(rocPixels):
            rpcs.push_back(vector[pixelConfig]())
            for pc in rocPixel:
                rpcs.at(idx).push_back(<pixelConfig> pc.thisptr[0])

        if rocI2C is not None:
            for i in rocI2C:
                i2c.push_back(i)
            return self.thisptr.initDUT(hubs, tbmtype, td, roctype,rd,rpcs,i2c)
        else:
            return self.thisptr.initDUT(hubs, tbmtype, td, roctype,rd,rpcs)

    def getVersion(self):
        return self.thisptr.getVersion()
    def testAllPixels(self, bool enable, rocid = None):
        if rocid is not None:
            self.thisptr._dut.testAllPixels(enable,rocid)
        else:
            self.thisptr._dut.testAllPixels(enable)

    def getTbmDACs(self, int tbmid):
        r = self.thisptr._dut.getTbmDACs(tbmid)
        return {tup.first: tup.second for tup in r}
  
    def getRocDACs(self, int rocid):
        r = self.thisptr._dut.getDACs(rocid)
        return {tup.first: tup.second for tup in r}
    def getDACs(self, int rocid):
        return self.getRocDACs(rocid)
  
    def updateTrimBits(self, trimming, int rocid):
        cdef vector[pixelConfig] v
        cdef pixelConfig pc
        #for idx, col, row, trim in enumerate(trimming):
        for line in xrange(len(trimming)):
            pc = pixelConfig(trimming[line][0][0], trimming[line][1][0], trimming[line][2][0])
            v.push_back(pc)
        self.thisptr._dut.updateTrimBits(v, rocid)

    def info(self):
        self.thisptr._dut.info()

    def setROCEnable(self, int rocid, bool enable):
        self.thisptr._dut.setROCEnable(rocid, enable)

    def setTBMEnable(self, int tbmid, bool enable):
        self.thisptr._dut.setTBMEnable(tbmid, enable)

    def testPixel(self, int col, int row, bool enable, rocid = None):
        if rocid is not None:
            self.thisptr._dut.testPixel(col, row, enable,rocid)
        else:
            self.thisptr._dut.testPixel(col, row, enable)
    def maskAllPixels(self, bool enable, rocid = None):
        if rocid is not None:
            self.thisptr._dut.maskAllPixels(enable,rocid)
        else:
            self.thisptr._dut.maskAllPixels(enable)
    def maskPixel(self, int col, int row, bool enable, rocid = None):
        if rocid is not None:
            self.thisptr._dut.maskPixel(col, row, enable,rocid)
        else:
            self.thisptr._dut.maskPixel(col, row, enable)

    def getNMaskedPixels(self, rocid = None):
        if rocid is not None:
            return self.thisptr._dut.getNMaskedPixels(rocid)
        else:
            return self.thisptr._dut.getNMaskedPixels()

    def getMaskedPixels(self, rocid = None):
        cdef vector[pixelConfig] rpcs
        if rocid is not None:
            rpcs = self.thisptr._dut.getMaskedPixels(rocid)
        else:
            rpcs = self.thisptr._dut.getMaskedPixels()
        pixelconfigs = list()
        for p in rpcs:
            pxc = PixelConfig()
            pxc.fill(p)
            pixelconfigs.append(pxc)
        return pixelconfigs

    def getNEnabledPixels(self, rocid = None):
        if rocid is not None:
            return self.thisptr._dut.getNEnabledPixels(rocid)
        else:
            return self.thisptr._dut.getNEnabledPixels()

    def getEnabledPixels(self, rocid = None):
        cdef vector[pixelConfig] rpcs
        if rocid is not None:
            rpcs = self.thisptr._dut.getEnabledPixels(rocid)
        else:
            rpcs = self.thisptr._dut.getEnabledPixels()
        pixelconfigs = list()
        for p in rpcs:
            pxc = PixelConfig()
            pxc.fill(p)
            pixelconfigs.append(pxc)
        return pixelconfigs

    def getNEnabledTbms(self):
        return self.thisptr._dut.getNEnabledTbms()
    def getNEnabledRocs(self):
        return self.thisptr._dut.getNEnabledRocs()

    def getEnabledRocI2Caddr(self):
        cdef vector[uint8_t] rpcs
        rpcs = self.thisptr._dut.getEnabledRocI2Caddr()
        roci2c = list()
        for r in rpcs:
            roci2c.append(r)
        return roci2c

    def getEnabledRocIDs(self):
        cdef vector[uint8_t] rpcs
        rpcs = self.thisptr._dut.getEnabledRocIDs()
        rocids = list()
        for r in rpcs:
            rocids.append(r)
        return rocids

    def getNTbms(self):
        return self.thisptr._dut.getNTbmCores()
    def getNRocs(self):
        return self.thisptr._dut.getNRocs()
    def getTbmType(self):
        return self.thisptr._dut.getTbmType()
    def getRocType(self):
        return self.thisptr._dut.getRocType()
    #def programDUT(self):
        #return self.thisptr.programDUT()
    def status(self):
        return self.thisptr.status()
    def flashTB(self, string filename):
        return self.thisptr.flashTB(filename)
    def getTBia(self):
        return float(self.thisptr.getTBia())
    def getTBva(self):
        return float(self.thisptr.getTBva())
    def getTBid(self):
        return float(self.thisptr.getTBid())
    def getTBvd(self):
        return float(self.thisptr.getTBvd())
    def HVoff(self):
        self.thisptr.HVoff()
    def HVon(self):
        self.thisptr.HVon()
    def Poff(self):
        self.thisptr.Poff()
    def Pon(self):
        self.thisptr.Pon()
    def SignalProbe(self, string probe, string name, int channel = 0):
        return self.thisptr.SignalProbe(probe, name, channel)
    def setDAC(self, string dacName, uint8_t dacValue, rocid = None):
        if rocid is None:
            return self.thisptr.setDAC(dacName, dacValue)
        else:
            return self.thisptr.setDAC(dacName, dacValue, rocid)
    def getDACRange(self, string dacName):
        return self.thisptr.getDACRange(dacName)
    def setTbmReg(self, string regName, uint8_t regValue, tbmid = None):
        if tbmid is None:
            return self.thisptr.setTbmReg(regName, regValue)
        else:
            return self.thisptr.setTbmReg(regName, regValue, tbmid)
    def getPulseheightVsDAC(self, string dacName, int dacStep, int dacMin, int dacMax, int flags = 0, int nTriggers = 16):
        cdef vector[pair[uint8_t, vector[pixel]]] r
        r = self.thisptr.getPulseheightVsDAC(dacName, dacStep, dacMin, dacMax, flags, nTriggers)
        dac_steps = list()
        for d in xrange(r.size()):
            pixels = list()
            for pix in range(r[d].second.size()):
                p = r[d].second[pix]
                px = Pixel()
                px.fill(p)
                pixels.append(px)
            dac_steps.append(pixels)
        return numpy.array(dac_steps)

    def getEfficiencyVsDAC(self, string dacName, int dacStep, int dacMin, int dacMax, int flags = 0, int nTriggers = 16):
        cdef vector[pair[uint8_t, vector[pixel]]] r
        r = self.thisptr.getEfficiencyVsDAC(dacName, dacStep, dacMin, dacMax, flags, nTriggers)
        dac_steps = list()
        for d in xrange(r.size()):
            pixels = list()
            for pix in range(r[d].second.size()):
                p = r[d].second[pix]
                px = Pixel()
                px.fill(p)
                pixels.append(px)
            dac_steps.append(pixels)
        return numpy.array(dac_steps)

    def getEfficiencyVsDACDAC(self, string dac1name, uint8_t dac1step, uint8_t dac1min, uint8_t dac1max, string dac2name, uint8_t dac2step, uint8_t dac2min, uint8_t dac2max, uint16_t flags = 0, uint32_t nTriggers=16):
        cdef vector[pair[uint8_t, pair[uint8_t, vector[pixel]]]] r
        r = self.thisptr.getEfficiencyVsDACDAC(dac1name, dac1step, dac1min, dac1max, dac2name, dac2step, dac2min, dac2max, flags, nTriggers)
        # Return the linearized matrix with all pixels:
        dac_steps = list()
        for d in xrange(r.size()):
            pixels = list()
            for pix in xrange(r[d].second.second.size()):
                p = r[d].second.second[pix]
                px = Pixel()
                px.fill(p)
                pixels.append(px)
            dac_steps.append(pixels)
        return numpy.array(dac_steps)

    def getThresholdVsDAC(self, string dac1Name, uint8_t dac1Step, uint8_t dac1Min, uint8_t dac1Max, string dac2Name, uint8_t dac2Step, uint8_t dac2Min, uint8_t dac2Max, threshold, uint16_t flags = 0, uint32_t nTriggers=16):
        cdef vector[pair[uint8_t, vector[pixel]]] r
        r = self.thisptr.getThresholdVsDAC(dac1Name, dac1Step, dac1Min, dac1Max, dac2Name, dac2Step, dac2Min, dac2Max, threshold, flags, nTriggers)
        dac_steps = list()
        for d in xrange(r.size()):
            pixels = list()
            for pix in range(r[d].second.size()):
                p = r[d].second[pix]
                px = Pixel()
                px.fill(p)
                pixels.append(px)
            dac_steps.append(pixels)
        return numpy.array(dac_steps)

    def getPulseheightVsDACDAC(self, string dac1name, uint8_t dac1step, uint8_t dac1min, uint8_t dac1max, string dac2name, uint8_t dac2step, uint8_t dac2min, uint8_t dac2max, uint16_t flags = 0, uint32_t nTriggers=16):
        cdef vector[pair[uint8_t, pair[uint8_t, vector[pixel]]]] r
        r = self.thisptr.getPulseheightVsDACDAC(dac1name, dac1step, dac1min, dac1max, dac2name, dac2step, dac2min, dac2max, flags, nTriggers)
        # Return the linearized matrix with all pixels:
        dac_steps = list()
        for d in xrange(r.size()):
            pixels = list()
            for pix in xrange(r[d].second.second.size()):
                p = r[d].second.second[pix]
                px = Pixel()
                px.fill(p)
                pixels.append(px)
            dac_steps.append(pixels)
        return numpy.array(dac_steps)

    def getPulseheightMap(self, int flags, int nTriggers):
        cdef vector[pixel] r
        r = self.thisptr.getPulseheightMap(flags, nTriggers)
        pixels = list()
        for p in r:
            px = Pixel()
            px.fill(p)
            pixels.append(px)
        return pixels

    def getEfficiencyMap(self, int flags, int nTriggers):
        cdef vector[pixel] r
        r = self.thisptr.getEfficiencyMap(flags, nTriggers)
        pixels = list()
        for p in r:
            px = Pixel()
            px.fill(p)
            pixels.append(px)
        return pixels

    def getThresholdMap(self, string dacName, uint8_t dacStep, uint8_t dacMin, uint8_t dacMax, uint8_t threshold, int flags, int nTriggers):
        cdef vector[pixel] r
        r = self.thisptr.getThresholdMap(dacName, dacStep, dacMin, dacMax, threshold, flags, nTriggers)
        pixels = list()
        for p in r:
            px = Pixel()
            px.fill(p)
            pixels.append(px)
        return pixels

    def setExternalClock(self, bool enable):
        return self.thisptr.setExternalClock(enable)

    def setClockStretch(self, uint8_t src, uint16_t delay, uint16_t width):
        self.thisptr.setClockStretch(src, delay, width)

    def setSignalMode(self, string signal, string mode, uint8_t speed):
        self.thisptr.setSignalMode(signal, mode, speed)

    def daqStart(self, flags = None):
        if flags is not None:
            return self.thisptr.daqStart(flags)
        else:
            return self.thisptr.daqStart(0)

    def daqStatus(self):
        return self.thisptr.daqStatus()

    def daqClear(self):
        self.thisptr.daqClear()

    def daqTriggerSource(self, string source, uint32_t period = 0):
        return self.thisptr.daqTriggerSource(source, period)

    def daqSingleSignal(self, string signal):
        return self.thisptr.daqSingleSignal(signal)

    def daqTrigger(self, uint32_t nTrig, uint16_t period = 0):
        self.thisptr.daqTrigger(nTrig,period)

    def daqTriggerLoop(self, uint16_t period):
        self.thisptr.daqTriggerLoop(period)

    def daqTriggerLoopHalt(self):
        self.thisptr.daqTriggerLoopHalt()

    def daqGetEvent(self):
        cdef Event r
        r = self.thisptr.daqGetEvent()
        p = PxEvent()
        p.clone(r)
        return p

    def daqGetEventBuffer(self):
        cdef vector[Event] r
        r = self.thisptr.daqGetEventBuffer()
        pixelevents = list()
        for event in r:
            p = PxEvent()
            p.clone(event)
            pixelevents.append(p)
        return pixelevents

    def daqGetRawEvent(self):
        cdef rawEvent r
        hits = []
        r = self.thisptr.daqGetRawEvent()
        for i in range(r.data.size()):
            hits.append(r.data[i])
        return hits

    def daqGetBuffer(self):
        cdef vector[uint16_t] r
        r = self.thisptr.daqGetBuffer()
        return r

    def daqGetRawEventBuffer(self):
        # Since we're just returning the 16bit ints as rawEvent in python,
        # this is the same as dqGetBuffer:
        return self.thisptr.daqGetBuffer()

    def daqGetReadback(self):
        cdef vector[vector[uint16_t]] r
        r = self.thisptr.daqGetReadback()
        return r

    def daqGetXORsum(self, uint8_t channel):
        cdef vector[uint8_t] r
        r = self.thisptr.daqGetXORsum(channel)
        return r

    def daqStop(self):
        return self.thisptr.daqStop()

    def getStatistics(self):
        cdef statistics r
        r = self.thisptr.getStatistics()
        s = Statistics()
        s.c_clone(r)
        return s

    def setReportingLevel(self, string logLevel):
        self.thisptr.setReportingLevel(logLevel)

    def getReportingLevel(self):
        return self.thisptr.getReportingLevel()

cimport regdict
cdef class PyRegisterDictionary:
    cdef regdict.RegisterDictionary *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = regdict.getInstance()
    def __dealloc__(self):
        self.thisptr = NULL
    def getAllROCNames(self):
        names = []
        cdef vector[string] v = self.thisptr.getAllROCNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        return names
    def getAllDTBNames(self):
        names = []
        cdef vector[string] v = self.thisptr.getAllDTBNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        return names
    def getAllTBMNames(self):
        names = []
        cdef vector[string] v = self.thisptr.getAllTBMNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        return names

cimport probedict
cdef class PyProbeDictionary:
    cdef probedict.ProbeDictionary *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = probedict.getInstance()
    def __dealloc__(self):
        self.thisptr = NULL
    def getAllAnalogNames(self):
        names = []
        cdef vector[string] v = self.thisptr.getAllAnalogNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        return names
    def getAllDigitalNames(self):
        names = []
        cdef vector[string] v = self.thisptr.getAllDigitalNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        return names
    def getAllNames(self):
        names = []
        cdef vector[string] v = self.thisptr.getAllDigitalNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        v = self.thisptr.getAllAnalogNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        return names

cimport triggerdict
cdef class PyTriggerDictionary:
    cdef triggerdict.TriggerDictionary *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = triggerdict.getInstance()
    def __dealloc__(self):
        self.thisptr = NULL
    def getAllNames(self):
        names = []
        cdef vector[string] v = self.thisptr.getAllNames()
        for i in xrange(v.size()):
            names.append(v.at(i))
        return names

