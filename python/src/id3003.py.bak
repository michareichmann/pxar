#!/usr/bin/python2.6
## @file
## Implements the specific x-ray generator Iso Debyeflex 3003 from
## General Electric
## @ingroup xrayClient

import sys
import struct
import serial
import time

debug = False

## Implementation of the General Electric Iso Debyeflex 3003 x-ray
# generator
#
# This class is an implementation of the x-ray generator Iso Debyeflex
# 3003 from General Electric. It uses a serial interface to communicate
# with the device.
# \ingroup xrayClient
class id3003_xray_generator():
	## Constructor that takes the file name of the serial
	## device as argument
	##
	## The constructor takes the serial device file name
	## and creates a serial handle which is not opened
	## @param device Serial device file name
	def __init__(self, device):
		## Serial device file name
		self.serialdevice = device
		## Serial device handle
		self.serial = serial.Serial(device)
	## Destructor for clean up
	##
	## The destructor closes the serial connection
	## to the id3003 x-ray generator. It does not
	## change the state of the device.
	def __del__(self):
		self.serial.close()
	def is_open(self):
		return self.serial.isOpen()
	def test_communication(self):
		answer = self.serial_command_with_response("ID\n")
		return answer[:23] == "*- ISO-DEBYEFLEX 3003 -"

	def set_voltage(self, kV):
		if type(kV) != int or kV < 2 or kV > 60:
			return False
		command = "SV:{0:02d}\n".format(int(kV))
		self.serial_command_no_response(command)
		return self.get_voltage() == kV
	def get_voltage(self):
		answer = self.serial_command_with_response("VN\n")
		return int(answer[6:8])
	def set_current(self, mA):
		if type(mA) != int or mA < 2 or mA > 80:
			return False
		command = "SC:{0:02d}\n".format(int(mA))
		self.serial_command_no_response(command)
		return self.get_current() == mA
	def get_current(self):
		answer = self.serial_command_with_response("CN\n")
		return int(answer[6:8])
	def set_hv(self, on):
		if on:
			self.serial_command_no_response("HV:1\n")
		else:
			self.serial_command_no_response("HV:0\n")

		if self.get_hv() != on:
			return False

		if not on:
			return True

		# When turning on wait for voltage and current to settle
		answer = self.serial_command_with_response("VN\n")
		nominal_voltage = int(answer[1:-1])
		answer = self.serial_command_with_response("CN\n")
		nominal_current = int(answer[1:-1])
		actual_voltage = 0
		actual_current = 0
		trials = 0
		while actual_voltage != nominal_voltage or actual_current != nominal_current:
			if trials > 100:
				self.set_hv(0)
				return False
			answer = self.serial_command_with_response("VA\n")
			actual_voltage = int(answer[1:-1])
			answer = self.serial_command_with_response("CA\n")
			actual_current = int(answer[1:-1])
			trials += 1
			time.sleep(0.2)

		return True
	def get_hv(self):
		answer = self.serial_command_with_response("SR:01\n")
		return bool(int(answer[1:-1]) & (1 << 6))
	def set_beam_shutter(self, beamno, on):
		if beamno < 1 or beamno > 4:
			return False
		if on:
			command = "OS:{0:d}\n".format(beamno)
		else:
			command = "CS:{0:d}\n".format(beamno)
		self.serial_command_no_response(command)
		print "Waiting for confirmation ...",
		sys.stdout.flush()
		time.sleep(2)
		answer = ""
		if beamno == 1 or beamno == 2:
			answer = self.serial_command_with_response("SR:03\n")
		elif beamno == 3 or beamno == 4:
			answer = self.serial_command_with_response("SR:04\n")
		bit = 0
		if beamno == 1 or beamno == 3:
			bit = 6
		elif beamno == 2 or beamno == 4:
			bit = 2
		if bool(int(answer[1:-1]) & (1 << bit)) == bool(on):
			print "OK"
		else:
			print "failed"
		return bool(int(answer[1:-1]) & (1 << bit)) == bool(on)

	## Sends a serial command without response
	##
	## Sends a serial command (an ASCII string) to the device
	## and does not wait or read a response
	## @param command The ASCII command to be sent
	## @return Does not return any value
	def serial_command_no_response(self, command):
		if debug:
			print command.strip()
		self.serial.write(unicode(command))

	## Sends a serial command and waits and returns the response
	##
	## Sends a serial command (an ASCII string) to the device
	## and waits for a response. The response is then returned
	## @param command The ASCII command to be sent
	## @return The serial response (ASCII format)
	def serial_command_with_response(self, command):
		if debug:
			print command.strip()
		self.serial.write(command)
		self.serial.flush()
		answer = ""
		while len(answer) == 0 or answer[-1] != '\r':
			answer += self.serial.read(1)
		if debug:
			print ">" + answer.strip() + "<"
		return answer

gen = id3003_xray_generator("/dev/ttyID3003")

if len(sys.argv) == 4 and (sys.argv[1] == "open" or sys.argv[1] == "close") and sys.argv[2] == "shutter":
	try:
		number = int(sys.argv[3])
	except:
		print "Error: Invalid shutter argument!"
		sys.exit(1)
	if number < 0 or number > 3:
		print "Error: Invalid shutter argument!"
	if sys.argv[1] == "open":
		print "Opening shutter", number
		gen.set_beam_shutter(number, True)
	else:
		print "Closing shutter", number
		gen.set_beam_shutter(number, False)

elif len(sys.argv) == 4 and sys.argv[1] == "set" and (sys.argv[2] == "voltage" or sys.argv[2] == "current"):
	try:
		value = int(sys.argv[3])
	except:
		print "Error: Invalid value!"
		sys.exit(1)
	if sys.argv[2] == "voltage":
		print "Setting voltage to", value, "kV"
		gen.set_voltage(value)
	else:
		print "Setting current to", value, "mA"
		gen.set_current(value)
elif len(sys.argv) == 3 and sys.argv[1] == "get" and (sys.argv[2] == "voltage" or sys.argv[2] == "current"):
	if sys.argv[2] == "voltage":
		print gen.get_voltage()
	else:
		print gen.get_current()
elif len(sys.argv) == 3 and sys.argv[1] == "get" and sys.argv[2] == "hv":
	print gen.get_hv()
elif len(sys.argv) == 4 and sys.argv[1] == "set" and sys.argv[2] == "hv" and (sys.argv[3] == "on" or sys.argv[3] == "off"):
	if sys.argv[3] == "on":
		gen.set_hv(True)
	else:
		gen.set_hv(False)
else:
	print "Error: Argument error!"
