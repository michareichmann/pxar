#include "datatypes.h"
#include "helper.h"
#include "log.h"
#include "exceptions.h"
#include "constants.h"

namespace pxar {


  void pixel::decodeRaw(uint32_t raw, bool invert) {
    // Get the pulse height:
    setValue(static_cast<double>((raw & 0x0f) + ((raw >> 1) & 0xf0)));
    if((raw & 0x10) > 0) {
      LOG(logDEBUGAPI) << "invalid pulse-height fill bit from raw value of "<< std::hex << raw << std::dec << ": " << *this;
      throw DataInvalidPulseheightError("Error decoding pixel raw value");
    }

    // Decode the pixel address
    int r2 =    (raw >> 15) & 7;
    if(invert) { r2 ^= 0x7; }
    int r1 = (raw >> 12) & 7;
    if(invert) { r1 ^= 0x7; }
    int r0 = (raw >>  9) & 7;
    if(invert) { r0 ^= 0x7; }
    int r = r2*36 + r1*6 + r0;
    _row = 80 - r/2;
    _column = 2*(((raw >> 21) & 7)*6 + ((raw >> 18) & 7)) + (r&1);

    // Perform range checks:
    if(_row >= ROC_NUMROWS || _column >= ROC_NUMCOLS) {
      LOG(logDEBUGAPI) << "Invalid pixel from raw value of "<< std::hex << raw << std::dec << ": " << *this;
      if(_row == ROC_NUMROWS) throw DataCorruptBufferError("Error decoding pixel raw value");
      else throw DataInvalidAddressError("Error decoding pixel raw value");
    }
  }

  void pixel::decodeLinear(uint32_t raw) {
    // Get the pulse height:
    setValue(static_cast<double>((raw & 0x0f) + ((raw >> 1) & 0xf0)));
    if((raw & 0x10) > 0) {
      LOG(logDEBUGAPI) << "invalid pulse-height fill bit from raw value of "<< std::hex << raw << std::dec << ": " << *this;
      throw DataInvalidPulseheightError("Error decoding pixel raw value");
    }

    // Perform checks on the fill bits:
    if((raw & 0x1000) > 0 || (raw & 0x100000) > 0) {
      LOG(logDEBUGAPI) << "invalid address fill bit from raw value of "<< std::hex << raw << std::dec << ": " << *this;
      throw DataInvalidAddressError("Error decoding pixel raw value");
    }

    // Decode the pixel address
    _column = ((raw >> 17) & 0x07) + ((raw >> 18) & 0x38);
    _row = ((raw >> 9) & 0x07) + ((raw >> 10) & 0x78);

    // Perform range checks:
    if(_row >= ROC_NUMROWS || _column >= ROC_NUMCOLS) {
      LOG(logDEBUGAPI) << "Invalid pixel from raw value of "<< std::hex << raw << std::dec << ": " << *this;
      if(_row == ROC_NUMROWS) throw DataCorruptBufferError("Error decoding pixel raw value");
      else throw DataInvalidAddressError("Error decoding pixel raw value");
    }
  }

  uint8_t pixel::translateLevel(int16_t x, int16_t level0, int16_t level1, int16_t levelS) {
    uint8_t level = uint8_t((x + level1 + levelS) / level1);
    return int(level) > 5 ? uint8_t(5) : level;
  }

  uint8_t pixel::translateLevel(int16_t x, std::vector<float> thresholds, uint8_t lastLevel, bool adjust) {

    int16_t level = adjust ? adjustLevel(expandSign(x), lastLevel, thresholds) : expandSign(x);
    for (uint8_t i(0); i < thresholds.size(); i++)
      if (level < thresholds.at(i))
        return i;
    return 5;
  }

  int16_t pixel::adjustLevel(int16_t analogue, uint8_t lastLevel, std::vector<float> thresholds) {
    float diff =  5 * (lastLevel - translateLevel(analogue, thresholds, lastLevel, false));
    return int16_t(analogue - diff);
  }

  void pixel::decodeAnalog(std::vector<uint16_t> analog, std::vector<float> thresholds) {

    /** Get the pulse height */
    setValue(static_cast<double>(expandSign(analog.back() & 0x0fff)));

    /** Get the column and row */
    uint8_t c1 = translateLevel(analog.at(0), thresholds, 5);
    uint8_t c0 = translateLevel(analog.at(1), thresholds, c1);
    int c  = c1*6 + c0;

    uint8_t r2 = translateLevel(analog.at(2), thresholds, c0);
    uint8_t r1 = translateLevel(analog.at(3), thresholds, r2);
    uint8_t r0 = translateLevel(analog.at(4), thresholds, r1);
    int r  = (r2*6 + r1)*6 + r0;

    _row = uint8_t(80 - r/2);
    _column = uint8_t(2*c + (r&1));

    /** Perform range checks */
    if(_row >= ROC_NUMROWS || _column >= ROC_NUMCOLS) {
      LOG(logDEBUGAPI) << "Invalid pixel from levels "<< listVector(analog) << ": " << *this;
      throw DataInvalidAddressError("Error decoding pixel address");
    }
  }

  void pixel::decodeAnalog(std::vector<uint16_t> analog, int16_t ultrablack, int16_t black) {
    // Check pixel data length:
    if(analog.size() != 6) {
      LOG(logDEBUGAPI) << "Received wrong number of data words for a pixel: " << analog.size();
      throw DataInvalidAddressError("Received wrong number of data words for a pixel");
    }

    // Calculate the levels:
    /** Changes by Micha*/
    int16_t level0 = black;
    int16_t level1 = (black - ultrablack)/4;
    int16_t levelS = level1/2;

    /** Get the pulse height */
    setValue(static_cast<double>(expandSign(analog.back() & 0x0fff) - level0));

    // Decode the pixel address
    int c1 = translateLevel(expandSign(analog.at(0)),level0,level1,levelS);
    int c0 = translateLevel(expandSign(analog.at(1)),level0,level1,levelS);
    int c  = c1*6 + c0;

    int r2 = translateLevel(expandSign(analog.at(2)),level0,level1,levelS);
    int r1 = translateLevel(expandSign(analog.at(3)),level0,level1,levelS);
    int r0 = translateLevel(expandSign(analog.at(4)) - 10,level0,level1,levelS);
    int r  = (r2*6 + r1)*6 + r0;

//    std::cout << ultrablack << ", " << black << ", ";
//    for (int i(0);i<5;i++)
//      std::cout << expandSign(analog.at(i)) << ", ";
//    std::cout << std::endl;

    _row = 80 - r/2;
    _column = 2*c + (r&1);

    /**Output by Micha*/
    std::stringstream ss;
    ss << "AnalogLevels: ";
    ss<<(int)_column<<" "<<(int)_row << "\t";
    for (unsigned i(0); i<analog.size(); i++)
        ss << analog[i] << " ";
    ss << "\t" << c1 <<" "<<c0<<" "<<r2<<" "<<r1<<" "<< r0;
    LOG(logDEBUGAPI)<<ss.str() << " ";

    /** Perform range checks:*/
    if(_row >= ROC_NUMROWS || _column >= ROC_NUMCOLS) {
      LOG(logDEBUGAPI) << "Invalid pixel from levels "<< listVector(analog) << ": " << *this;
      throw DataInvalidAddressError("Error decoding pixel address");
    }
  }

  uint32_t pixel::encode() {
    uint32_t raw = 0;
    // Set the pulse height:
    raw = ((static_cast<int>(value()) & 0xf0) << 1) + (static_cast<int>(value()) & 0xf);

    // Encode the pixel address
    int r = 2*(80 - _row);
    raw |= ((r/36) << 15);
    raw |= (((r%36)/6) << 12);
    raw |= (((r%36)%6 + _column%2) << 9);

    int dcol = _column/2;
    raw |= ((dcol)/6 << 21);
    raw |= (((dcol%6)) << 18);

    LOG(logDEBUGPIPES) << "Pix  " << static_cast<int>(_column) << "|"
		       << static_cast<int>(_row) << " = "
		       << dcol << "/" << r << " = "
		       << dcol/6 << " " << dcol%6 << " "
		       << r/36 << " " << (r%36)/6 << " " << (r%36)%6;

    // Return the 24 bits belonging to the pixel:
    return (raw & 0x00ffffff);
  }

  void Event::printHeader() {
    LOG(logINFO) << "Header content: 0x" << std::hex << header << std::dec;
    LOG(logINFO) << "\t Event ID \t" << static_cast<int>(this->triggerCount());
    LOG(logINFO) << "\t Data ID " << static_cast<int>(this->dataID()) 
		       << " Value " << static_cast<int>(this->dataValue());
  }

  void Event::printTrailer() {
    LOG(logINFO) << "Trailer content: 0x" << std::hex << trailer << std::dec;
    LOG(logINFO) << "\t Token Pass \t" << textBool(this->hasTokenPass());
    LOG(logINFO) << "\t Reset TBM \t" << textBool(this->hasResetTBM());
    LOG(logINFO) << "\t Reset ROC \t" << textBool(this->hasResetROC());
    LOG(logINFO) << "\t Sync Err \t" << textBool(this->hasSyncError());
    LOG(logINFO) << "\t Sync Trigger \t" << textBool(this->hasSyncTrigger());
    LOG(logINFO) << "\t ClearTrig Cnt \t" << textBool(this->hasClearTriggerCount());
    LOG(logINFO) << "\t Cal Trigger \t" << textBool(this->hasCalTrigger());
    LOG(logINFO) << "\t Stack Full \t" << textBool(this->stackFull());

    LOG(logINFO) << "\t Auto Reset \t" << textBool(this->hasAutoReset());
    LOG(logINFO) << "\t PKAM Reset \t" << textBool(this->hasPkamReset());
    LOG(logINFO) << "\t Stack Count \t" << static_cast<int>(this->stackCount());
  }

  void statistics::dump() {
    // Print out the full statistics:
    LOG(logINFO) << "Decoding statistics:";
    LOG(logINFO) << "  General information:";
    LOG(logINFO) << "\t 16bit words read:         " << this->info_words_read();
    LOG(logINFO) << "\t valid events total:       " << this->info_events_total();
    LOG(logINFO) << "\t empty events:             " << this->info_events_empty();
    LOG(logINFO) << "\t valid events with pixels: " << this->info_events_valid();
    LOG(logINFO) << "\t valid pixel hits:         " << this->info_pixels_valid();
    LOG(logINFO) << "  Event errors: \t           " << this->errors_event();
    LOG(logINFO) << "\t start marker:             " << this->errors_event_start();
    LOG(logINFO) << "\t stop marker:              " << this->errors_event_stop();
    LOG(logINFO) << "\t overflow:                 " << this->errors_event_overflow();
    LOG(logINFO) << "\t invalid 5bit words:       " << this->errors_event_invalid_words();
    LOG(logINFO) << "\t invalid XOR eye diagram:  " << this->errors_event_invalid_xor();
    LOG(logINFO) << "  TBM errors: \t\t           " << this->errors_tbm();
    LOG(logINFO) << "\t flawed TBM headers:       " << this->errors_tbm_header();
    LOG(logINFO) << "\t flawed TBM trailers:      " << this->errors_tbm_trailer();
    LOG(logINFO) << "\t event ID mismatches:      " << this->errors_tbm_eventid_mismatch();
    LOG(logINFO) << "  ROC errors: \t\t           " << this->errors_roc();
    LOG(logINFO) << "\t missing ROC header(s):    " << this->errors_roc_missing();
    LOG(logINFO) << "\t misplaced readback start: " << this->errors_roc_readback();
    LOG(logINFO) << "  Pixel decoding errors:\t   " << this->errors_pixel();
    LOG(logINFO) << "\t pixel data incomplete:    " << this->errors_pixel_incomplete();
    LOG(logINFO) << "\t pixel address:            " << this->errors_pixel_address();
    LOG(logINFO) << "\t pulse height fill bit:    " << this->errors_pixel_pulseheight();
    LOG(logINFO) << "\t buffer corruption:        " << this->errors_pixel_buffer_corrupt();
  }

  void statistics::clear() {
    m_info_words_read = 0;
    m_info_events_empty = 0;
    m_info_events_valid = 0;
    m_info_pixels_valid = 0;

    m_errors_event_start = 0;
    m_errors_event_stop = 0;
    m_errors_event_overflow = 0;
    m_errors_event_invalid_words = 0;
    m_errors_event_invalid_xor = 0;

    m_errors_tbm_header = 0;
    m_errors_tbm_trailer = 0;
    m_errors_tbm_eventid_mismatch = 0;

    m_errors_roc_missing = 0;
    m_errors_roc_readback = 0;

    m_errors_pixel_incomplete = 0;
    m_errors_pixel_address = 0;
    m_errors_pixel_pulseheight = 0;
    m_errors_pixel_buffer_corrupt = 0;
  }

  statistics& operator+=(statistics &lhs, const statistics &rhs) {
    // Informational bits:
    lhs.m_info_words_read += rhs.m_info_words_read;
    lhs.m_info_events_empty += rhs.m_info_events_empty;
    lhs.m_info_events_valid += rhs.m_info_events_valid;
    lhs.m_info_pixels_valid += rhs.m_info_pixels_valid;

    // Event errors:
    lhs.m_errors_event_start += rhs.m_errors_event_start;
    lhs.m_errors_event_stop += rhs.m_errors_event_stop;
    lhs.m_errors_event_overflow += rhs.m_errors_event_overflow;
    lhs.m_errors_event_invalid_words += rhs.m_errors_event_invalid_words;
    lhs.m_errors_event_invalid_xor += rhs.m_errors_event_invalid_xor;

    // TBM errors:
    lhs.m_errors_tbm_header += rhs.m_errors_tbm_header;
    lhs.m_errors_tbm_trailer += rhs.m_errors_tbm_trailer;
    lhs.m_errors_tbm_eventid_mismatch += rhs.m_errors_tbm_eventid_mismatch;

    // ROC errors:
    lhs.m_errors_roc_missing += rhs.m_errors_roc_missing;
    lhs.m_errors_roc_readback += rhs.m_errors_roc_readback;

    // Pixel decoding errors:
    lhs.m_errors_pixel_incomplete += rhs.m_errors_pixel_incomplete;
    lhs.m_errors_pixel_address += rhs.m_errors_pixel_address;
    lhs.m_errors_pixel_pulseheight += rhs.m_errors_pixel_pulseheight;
    lhs.m_errors_pixel_buffer_corrupt += rhs.m_errors_pixel_buffer_corrupt;

    return lhs;
  }

} // namespace pxar
