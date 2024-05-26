#!/usr/bin/env python

from curses import wrapper
import curses, sys, time
import curses.ascii
import math, decimal
from enum import Enum
import numpy as np
import re

Dec = decimal.Decimal
gl_pi = Dec('3.1415926535897932384626433833')
gl_e = Dec('2.7182818284590452353602874714')
gl_prec = 3

def remove_exponent(d):
	"""Remove the exponent of a given Decimal."""
	decimal.getcontext().prec -= 3
	d = +d
	decimal.getcontext().prec += 3
	return d.quantize(Dec(1)) if d == d.to_integral() else d.normalize()

def numberform(number):
	"""Return a string with a representation of a given number."""
	def convert(n, exp):
		mantissa = str(remove_exponent(n * Dec('1e' + str(-exp))))
		exponent = ' * 10^' + str(exp)
		return mantissa + exponent
	def compare(n, exp):
		return Dec('0.5e'+str(exp)) < abs(n) <= Dec('500e'+str(exp))
	if compare(number, 0) or compare(number, -1) or number == 0:
		ans = str(remove_exponent(number))
	else:
		for l in (3, 6, 9, -3, -6, -9):
			if compare(number, l):
				ans = convert(number, l)
				break
		else:
			l = math.floor(Dec.log10(abs(number)))
			ans = convert(number, l)
	return ans

def mass(cp):
	def f(st):
		if st.isalpha():
			if st not in table:
				raise ValueError('incorrect compound name')
			return table[st]
		else:
			if len(st) < 2:
				raise ValueError('incorrect compound name')
			elif st[-2].isdigit():
				st, n = st[:-2], st[-2:]
			else:
				st, n = st[:-1], st[-1:]
			if st not in table:
				raise ValueError('incorrect compound name')
			return table[st] * int(n)
	
	table = {'H': 1, 'Li': 7, 'Na': 23, 'K': 39, 'Rb': 85, 'Cs': 133, 'Fr': 223, 
		'Be': 9, 'Mg': 24, 'Ca': 40, 'Sr': 88, 'Ba': 137, 'Ra': 226, 
		'Sc': 45, 'Y': 89, 'La': 139, 'Ac': 227, 
		'Ti': 48, 'Zr': 91, 'Hf': 178, 'Rf': 265, 
		'V': 51, 'Nb': 93, 'Ta': 181, 'Db': 268, 
		'Cr': 52, 'Mo': 96, 'W': 184, 'Sg': 271, 
		'Mn': 55, 'Tc': 98, 'Re': 186, 'Bh': 272, 
		'Fe': 56, 'Ru': 101, 'Os': 190, 'Hs': 270, 
		'Co': 59, 'Rh': 103, 'Ir': 192, 'Mt': 276, 
		'Ni': 59, 'Pd': 106, 'Pt': 195, 'Ds': 281, 
		'Cu': 64, 'Ag': 108, 'Au': 197, 'Rg': 280, 
		'Zn': 65, 'Cd': 112, 'Hg': 201, 'Cn': 285, 
		'B': 11, 'Al': 27, 'Ga': 70, 'In': 115, 'Tl': 204, 'Nh': 284, 
		'C': 12, 'Si': 28, 'Ge': 73, 'Sn': 119, 'Pb': 207, 'Fl': 289, 
		'N': 14, 'P': 31, 'As': 75, 'Sb': 122, 'Bi': 209, 'Mc': 288, 
		'O': 16, 'S': 32, 'Se': 79, 'Te': 128, 'Po': 209, 'Lv': 293, 
		'F': 19, 'Cl': 35.5, 'Br': 80, 'I': 127, 'At': 210, 'Ts': 294, 
		'He': 4, 'Ne': 20, 'Ar': 40, 'Kr': 84, 'Xe': 131, 'Rn': 222, 'Og': 294,
		'Ce': 140, 'Pr': 141, 'Nd': 144, 'Pm': 145, 'Sm': 150, 'Eu': 152, 'Gd': 157, 
		'Tb': 159, 'Dy': 162, 'Ho': 165, 'Er': 167, 'Tm': 169, 'Yb': 173, 'Lu': 175, 
		'Th': 232, 'Pa': 231, 'U': 238, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247, 
		'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262}
	
	cp = re.sub(r'([A-Z()\[\]*])', r' \1', cp)[1:]
	ls = cp.split()
	m = [0, 0, 0, 0]
	level = 1
	molec_num = 1
	for el in ls:
		if el in ('(', '['):
			level += 1
		elif ')' in el or ']' in el:
			n = 1
			if len(el) > 1:
				n = int(el[1:])
			m[level - 1] += m[level] * n
			m[level] = 0
			level -= 1
		elif el == '*':
			m[0] += m[1]
			m[1] = 0
		else:
			m[level] += f(el)
	m[0] += m[1]
	m[1] = 0
	if m[1] + m[2] + m[3] != 0:
		raise ValueError('incorrect compound name')
	return Dec(m[0])

def prime_fact(n):
	"""Return the prime factorisation of n."""
	primes = []
	ans = []
	n = int(n)
	x = 1
	power = 0
	while n % 2 == 0:
		n //= 2
		power += 1
	if power > 0:
		ans.append((2, power))
	while n > 1:
		x += 2
		power = 0
		is_prime = True
		if x > 10000000:
			raise ValueError('large number pf')
		for a in primes:
			if a * a > x:
				is_prime = True
				break
			if x % a == 0:
				is_prime = False
				break
		if not is_prime:
			continue
		primes.append(x)
		while n % x == 0:
			n //= x
			power += 1
		if power != 0:
			ans.append((x, power))
		if x * x > n and n != 1:
			ans.append((n, 1))
			return ans
	return ans

def prime_fact_show(n):
	toprint = ""
	ans = prime_fact(n)
	for el in ans[:-1]:
		if el[1] == 1:
			toprint += str(el[0]) + ' * '
		else:
			toprint += str(el[0]) + '^' + str(el[1]) + ' * '
	el = ans[-1]
	if el[1] == 1:
		toprint += str(el[0])
	else:
		toprint += str(el[0]) + '^' + str(el[1])
	return toprint

def fctr(n):
	ans = 1
	n = int(n)
	for t in range(n):
		ans *= (t + 1)
	return Dec(ans)

def P(n, k):
	ans = 1
	k = int(k)
	for a in range(k):
		ans *= (n - a)
	return Dec(ans)

def C(n, k):
	return P(n, k) / fctr(k)

class AttributeSet:
	"""A class for keeping attributes."""
	pass


class Unit:
	"""The creation of the Unit object and the related functionality."""
	
	standart_units = ['m', 's', 'kg']
	
	class InvalidOperationError(ArithmeticError):
		pass
	
	def __init__(self, value, units):
		"""The initialiser for the class.
		
		Arguments:
		value -- a number representing the numeric value of a quantity
		units -- a dictionary that matches units and their powers
		"""
		self.units = units
		self.value = value
		
	def usr_init(*kwards):
		"""Get a Unit without using a dictionary.
		
		>>> Unit.usr_init(5, 'm')
		Unit(5 m)
		"""
		ansdict = {}
		ansnum = Dec(1)
		for u in kwards:
			upow = 1
			if type(u) in (int, Dec):
				ansnum *= u
			else:
				if u[-1].isdigit():
					u, upow = u[:-1], int(u[-1])
				ansdict |= {u: upow}
		if ansdict == {}:
			return ansnum
		return Unit(ansnum, ansdict)
	
	def __str__(self):
		ans = ''
		for u in list(self.units):
			power = self.units[u]
			ans += str(u)
			if power != 1:
				ans += '^' + str(power)
			ans += '*'
		return ans[:-1]
	
	def smart_str(self, mode=0):
		"""Transform a unit into a string according to mode."""
		if min(self.units[u] for u in self.units) < 0:
			return self.fract_str(mode)
		elif mode == 2:
			return str(self.value) + ' ' + str(self)
		elif mode == 1:
			return numberform(self.value) + ' ' + str(self)
		elif self.units == {'s': 1}:
			return self.time_str()
		elif self.units == {'kg': 1}:
			return self.mass_str()
		elif self.units == {'m': 1}:
			return self.len_str()
		else:
			return numberform(self.value) + ' ' + str(self)
	
	def convert(self, n, exp):
		"""Calculate n * 10^exp."""
		return remove_exponent(n * Dec('1e' + str(-exp)))
	
	def compare(self, n, exp):
		"""Check if n is between 0.5 * 10^exp and 500 * 10^exp."""
		return Dec('0.5e'+str(exp)) < abs(n) <= Dec('500e'+str(exp))
	
	def find_exponent(self, value, units):
		"""Transform value to the first close unit in units.
		
		Arguments:
		value -- the value to be transformed
		units -- a tuple of tuples (exponent, name)
		def_unit -- the unit in which value is given.
		"""
		for l in units:
			n, s = l
			if self.compare(value, n):
				return str(self.convert(value, n)) + ' ' + s
		else:
			l = math.floor(math.log10(abs(value)))
			return str(self.convert(value, l)) + ' ' + units[0][1]
	
	def fract_str(self, mode):
		"""Represent the unit as a fraction."""
		num_ls = []
		denom_ls = []
		for unit in self.units:
			if self.units[unit] < 0:
				denom_ls.append(unit)
			else:
				num_ls.append(unit)
		if mode == 2:
			main = str(self.value)
		else:
			main = numberform(self.value)
		num_unit = Unit(1, {u: self.units[u] for u in num_ls})
		denom_unit = Unit(1, {u: -self.units[u] for u in denom_ls})
		num = str(num_unit)
		denom = str(denom_unit)
		if num == '':
			num = '1'
		l = (max(len(num), len(denom)) + 2)
		num = ' ' * ((l - len(num)) // 2 + len(main)) + num
		num += ' ' * ((l - len(num) + 1) // 2)
		denom = ' ' * ((l - len(denom)) // 2 + len(main)) + denom
		denom += ' ' * ((l - len(denom) + 1) // 2)
		main += '╶' + '─' * (l - 2) + '╴'
		return (num, main, denom)
	
	def time_str(self):
		"""Transform time from seconds into a practical measuring system."""
		ans = ''
		# time is measured in milliseconds
		v = int(self.value * 1000)
		year = 31557600000  # 365.25 days
		day, hour = 86400000, 3600000
		minute, second = 60000, 1000
		name_of = {year: 'year', day: 'day', hour: 'h'}
		name_of |= {minute: 'min', second: 's'}
		
		for period in (year, day, hour, minute, second):
			if v >= period:
				a = v // period
				if period < day or a in (1, -1):
					ans += str(a) + ' ' + name_of[period] + ' '
				else:
					ans += str(a) + ' ' + name_of[period] + 's '
				v %= period
		if v != 0:
			ans += str(v) + ' ms '
		
		return ans[:-1]
	
	def mass_str(self):
		"""Transform mass into a practical measuring system."""
		ans = ''
		v = self.value
		units = ((0, 'kg'), (3, 't'), (6, 'kt'), (9, 'Mt'), (-3, 'g'), (-6, 'mg'), (-9, 'μg'))
		if v == 0:
			return '0'
		elif Dec('0.01') <= v < 1:
			return str(self.convert(v, -3)) + ' g'
		else:
			return self.find_exponent(v, units)
	
	def len_str(self):
		"""Transform length into a practical measuring system."""
		ans = ''
		v = self.value
		units = ((0, 'm'), (3, 'km'), (-3, 'mm'), (-6, 'μm'), (-9, 'nm'))
		if v == 0:
			return '0'
		elif Dec('0.01') <= v < 1:
			return str(self.convert(v, -2)) + ' cm'
		elif v > 500000:
			v /= Dec(1000)
			l = math.floor(math.log10(abs(v)))
			return str(self.convert(v, l)) + f' * 10^{l} km'
		else:
			return self.find_exponent(v, units)
	
	def getpow(self, u):
		"""Return the power in which unit u is present in the object."""
		try:
			return self.units[u]
		except KeyError:
			return 0
	
	def __mul__(self, other):
		if type(other) == type(self):
			su = list(self.units | other.units)
			ansunit = {m: self.getpow(m) + other.getpow(m) for m in su}
			for m in list(ansunit):
				if ansunit[m] == 0:
					del ansunit[m]
			if ansunit == {}:
				return self.value * other.value
			return type(self)(self.value * other.value, ansunit)
		else:
			if isinstance(other, Unit):
				oth = other.value
			else:
				oth = other
			return type(self)(self.value * oth, self.units)
	
	def __truediv__(self, other):
		if type(other) == type(self):
			su = list(self.units | other.units)
			ansunit = {m: self.getpow(m) - other.getpow(m) for m in su}
			for m in list(ansunit):
				if ansunit[m] == 0:
					del ansunit[m]
			if ansunit == {}:
				return self.value / other.value
			return type(self)(self.value / other.value, ansunit)
		else:
			if isinstance(other, Unit):
				oth = other.value
			else:
				oth = other
			return type(self)(self.value / oth, self.units)
		
	def __pos__(self):
		return type(self)(self.value, self.units)
	
	def __add__(self, other):
		if type(other) is type(self):
			if self.units == other.units:
				return type(self)(self.value + other.value, self.units)
		raise Unit.InvalidOperationError('cannot add different units')
	
	def __neg__(self):
		return type(self)(-self.value, self.units)
	
	def __sub__(self, other):
		if type(other) is type(self):
			if self.units == other.units:
				return type(self)(self.value - other.value, self.units)
		raise Unit.InvalidOperationError('cannot substract different units')
	
	def __rmul__(self, other):
		return self * other
	
	def __rtruediv__(self, other):
		return other * Unit(1, {}) / self
	
	def __radd__(self, other):
		raise TypeError('cannot add different units')
	
	def __rsub__(self, other):
		raise TypeError('cannot add different units')
	
	def __round__(self, num):
		return type(self)(round(self.value, num), self.units)
	
	def __pow__(self, other, opt=None):
		if type(other) in [int, Dec]:
			ansunit = {a: self.units[a] * other for a in list(self.units)}
			return type(self)(self.value ** other, ansunit)
		else:
			raise Unit.InvalidOperationError('cannot raise to a (unit) power')
	
	def __repr__(self):
		if self.value == 1:
			return f'Unit({str(self)})'
		return f'Unit({str(self.value)}, {str(self)})'
	
	def si_units():
		"""Return a dictionary with units of SI."""
		
		# ways to initialise units
		one = lambda a: Unit.usr_init(Dec(1), a)
		derived = lambda a: Unit(Dec(1), a)
		U = Unit.usr_init
		
		unitdict = dict()  # store defined units here
		prefixlist = ['', 'da', 'h', 'k', '$', '$', 'M', '$', '$', 'G']
		prefixlist += ['n', '$', '$', 'μ', '$', '$', 'm', 'c', 'd']
		
		def add_units(u, name, explist=[9, 6, 3, 0, -3, -6, -9]):
			nonlocal unitdict
			for exp in explist:
				u2 = u * Dec('1e' + str(exp))
				name2 = prefixlist[exp] + name
				unitdict |= {name2: u2}
		
		gramm = U(Dec('1e-3'), 'kg')
		tonne = U(Dec('1e3'), 'kg')
		
		add_units(gramm, 'g', [3, 0, -3, -6, -9])
		add_units(tonne, 't', [9, 6, 3, 0])
		add_units(one('m'), 'm', [3, 0, -1, -2, -3, -6, -9])
		add_units(one('s'), 's', [0, -3, -6, -9])
		add_units(one('A'), 'A')
		add_units(one('K'), 'K')
		add_units(one('mol'), 'mol')
		
		litre = Unit(Dec('1e-3'), {'m': 3})
		herz = derived({'s': -1})
		newton = derived({'kg': 1, 'm': 1, 's': -2})
		pascal = derived({'kg': 1, 'm': -1, 's': -2})
		joule = derived({'kg': 1, 'm': 2, 's': -2})
		watt = derived({'kg': 1, 'm': 2, 's': -3})
		coulomb = derived({'s': 1, 'A': 1})
		volt = derived({'kg': 1, 'm': 2, 's': -3, 'A': -1})
		ohm = derived({'kg': 1, 'm': 2, 's': -3, 'A': -2})
		becquerrel = derived({'s': -1})
		gray = derived({'m': 2, 's': -1})
		
		add_units(litre, 'l', [0, -3])
		add_units(herz, 'Hz')
		add_units(newton, 'N')
		add_units(pascal, 'Pa', [9, 6, 3, 2, 0, -3, -6, -9])
		add_units(joule, 'J')
		add_units(watt, 'W')
		add_units(coulomb, 'C')
		add_units(volt, 'V')
		add_units(ohm, 'Ω')
		add_units(becquerrel, 'Bq')
		add_units(gray, 'Gy')
		
		unitdict |= {'min': U(Dec(60), 's'), 'h': U(Dec(3600), 's')}
		
		return unitdict


class Angle(Unit):
	"""The creation of the Angle object and the related functionality."""
	
	def usr_init(value, deg=False):
		"""Get an Angle without using a dictionary.
		
		>>> Angle.usr_init(5)
		Angle(5)
		>>> Angle.usr_init(1, 'deg')
		Angle(0.0175)
		"""
		if deg:
			return Angle(value / Dec(180) * gl_pi, {'ang': 1})
		else:
			return Angle(value, {'ang': 1})
		
	def init(self):
		"""Change the angle if it is bigger than 2pi to equivalent."""
		v = self.value
		while v > 6 and v >= gl_pi * 2:
			v -= self.pi * 2
		self.value = v
		
	def degree(self):
		"""Return the value of the angle in degrees."""
		return round(self.value * Dec(180) / gl_pi, gl_prec)
	
	def pirad(self):
		"""Return the ratio between the angle and pi."""
		return round(gl_pi / self.value, gl_prec)
	
	def __str__(self):
		return str(self.value)
	
	def __repr__(self):
		return 'Angle(' + str(round(self.value, 4)) + ')'
	
	def get_value(self):
		"""Return the value of the angle in radians."""
		if type(self) is Angle:
			return self.value
		else:
			return Dec(self)
	
	def cos(self):
		"""Return the cosine of the angle."""
		if type(self) is Dec:
			self = Angle.usr_init(self)
		x = Angle.get_value(self)
		decimal.getcontext().prec += 2
		i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
		while s != lasts:
			lasts = s
			i += 2
			fact *= i * (i-1)
			num *= x * x
			sign *= -1
			s += num / fact * sign
		decimal.getcontext().prec -= 2
		return +s
	
	def sin(self):
		"""Return the sine of the angle."""
		if type(self) is Dec:
			self = Angle.usr_init(self)
		x = Angle.get_value(self)
		decimal.getcontext().prec += 2
		i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
		while s != lasts:
			lasts = s
			i += 2
			fact *= i * (i-1)
			num *= x * x
			sign *= -1
			s += num / fact * sign
		decimal.getcontext().prec -= 2
		return +s
	
	def tg(self):
		"""Return the tangene of the angle."""
		if type(self) is Dec:
			self = Angle.usr_init(self)
		return self.sin() / self.cos()
	
	def ctg(self):
		"""Return the cotangene of the angle."""
		if type(self) is Dec:
			self = Angle.usr_init(self)
		return self.cos() / self.sin()
	
	def arcsin(x):
		"""Return an angle with given sine."""
		return Angle(Dec(math.asin(x)), {'ang': 1})
	
	def arccos(x):
		"""Return an angle with given cosine."""
		return Angle(Dec(math.acos(x)), {'ang': 1})
	
	def arctg(x):
		"""Return an angle with given tangene."""
		return Angle(Dec(math.atan(x)), {'ang': 1})


class Operator(Enum):
	"""Operator objects and the related functionality."""
	
	# NAME = (calc, pref, arg_num, left)
	PLUS = (lambda a, b: a + b, 0, 2, 0)
	MINUS = (lambda a, b: a - b, 0, 2, 0)
	MUL = (lambda a, b: a.dot(b) if type(a) is np.ndarray and type(b) is np.ndarray else a * b, 1, 2, 0)
	DIV = (lambda a, b: a / b, 1, 2, 0)
	SPACE = (lambda a, b: a * b, 2, 2, 1)
	SQRT = (lambda x: x ** Dec('0.5'), 3, 1, 1)
	CBRT = (lambda x: x ** (Dec('1') / Dec('3')), 3, 1, 1)
	POW = (lambda a, b: sum(a ** b) if isinstance(a, np.ndarray) else a ** b, 3, 2, 1)
	NEG = (lambda x: x.__neg__(), 3, 1, 1)
	FRBAR = (lambda a, b: a / b, 4, 2, 0)
	
	SEMICOLON = (lambda a: a, -1, 1, 1)
	
	SIN = (Angle.sin, 2, 1, 1)
	COS = (Angle.cos, 2, 1, 1)
	TG = (Angle.tg, 2, 1, 1)
	CTG = (Angle.ctg, 2, 1, 1)
	ARCSIN = (Angle.arcsin, 2, 1, 1)
	ARCCOS = (Angle.arccos, 2, 1, 1)
	ARCTG = (Angle.arctg, 2, 1, 1)
	ABS = (abs, 2, 1, 1)
	ANGLE = (Angle.usr_init, 2, 1, 1)
	LN = (Dec.ln, 2, 1, 1)
	PF = (prime_fact_show, 2, 1, 1)
	
	SINPOW = (lambda a, b: Angle.sin(b) ** a, 2, 2, 1)
	COSPOW = (lambda a, b: Angle.cos(b) ** a, 2, 2, 1)
	TGPOW = (lambda a, b: Angle.tg(b) ** a, 2, 2, 1)
	CTGPOW = (lambda a, b: Angle.ctg(b) ** a, 2, 2, 1)
	
	MASS = (mass, 5, 1, 1)
	UNIT = (Unit.usr_init, 5, 1, 1)
	
	FCTR = (fctr, 5, 1, 1)
	C = (lambda g: C(g[0], g[1]), 5, 1, 1)
	P = (lambda g: P(g[0], g[1]), 5, 1, 1)
	
	SUM = (sum, 5, 1, 1)
	AVERAGE = (lambda g: 0 if len(g) == 0 else sum(g) / Dec(len(g)), 5, 1, 1)
	NORM = (lambda g: sum(g**2) ** Dec('0.5'), 5, 1, 1)
	
	def __init__(self, calc, pref, arg_num, left):
		"""The initialiser for the class.
		
		Arguments:
		calc -- the function of the operator
		pref -- a number 0-4 representing the operators preference
		unary -- if truthy, the operator is unary
		left -- if truthy, a list of operators is evaluated right to left
		"""
		self.calc = calc
		self.pref = pref
		self.arg_num = arg_num
		self.left = bool(left)
	
	def __repr__(self):
		return self.name


class Calculator:
	"""The creation of the Calculator object and the related functionality."""
	
	operdict = {'+': Operator.PLUS,  '-': Operator.MINUS,
				'*': Operator.MUL,   ':': Operator.DIV,
				'^': Operator.POW,   '~': Operator.NEG,
				'√': Operator.SQRT,  '∛': Operator.CBRT,
				'/': Operator.FRBAR, '@': Operator.SPACE,
				';': Operator.SEMICOLON,
				'sin': Operator.SIN, 'arcsin': Operator.ARCSIN,
				'cos': Operator.COS, 'arccos': Operator.ARCCOS,
				'tg': Operator.TG,   'arctg': Operator.ARCTG,
				'ctg': Operator.CTG, 'cot': Operator.CTG,
				'tan': Operator.TG,  'arctan': Operator.ARCTG,
				'abs': Operator.ABS,
				'sin_': Operator.SINPOW, 'cos_': Operator.COSPOW,
				'tg_': Operator.TGPOW,   'ctg_': Operator.CTGPOW,
				'tan_': Operator.TGPOW,  'cot_': Operator.CTGPOW,
				'ln': Operator.LN,    'mass': Operator.MASS,
				'comb': Operator.C,   'perm': Operator.P, 
				'fctr': Operator.FCTR,'pf': Operator.PF,
				'unit': Operator.UNIT,'angle': Operator.ANGLE,
				'sum': Operator.SUM,  'average': Operator.AVERAGE,
				'norm': Operator.NORM}
	oper_regexp = r'+\-*:^~√∛/@$'
	helptext = [" " * 66,
	"# Welcome to the fractum manual!                                  ",
	"Fractum is a command line scientific calculator. It is designed to",
	"display mathematical equations as if they were written on paper.  ",
	"To do this, it has horizontal fraction bars and allows you to omit",
	"parentheses  when  calling  a  function  (more  on  this  below). ",
	"                                                                  ",
	"The manual is made up of several topics, marked with '#'.         ",
	"Some information is wrapped into columns for better readability.  ",
	"The definition of the column is on the top, marked with curly     ",
	"brackets. To isolate examples, parentheses are used (example).    ",
	"                                                                  ",
	"# Screen.                                                         ",
	"The calculator's screen is made up of five parts:                 ",
	"{part}                     {description}                          ",
	"top line                   name and version                       ",
	"answer status              type and form of the answer            ",
	"status line                key suggestions                        ",
	"main calculator area       equation and answer / manual           ",
	"three bottom lines         keyboard shortcuts                     ",
	"                                                                  ",
	"# Working with data                                               ",
	"Fractum supports 5 main data storing objects: numbers, vectors,   ",
	"units, angles and strings; and 2 main data transforming tokens:   ",
	"operators and functions. You will learn more about them in the    ",
	"following topics.                                                 ",
	"                                                                  ",
	"# Numbers                                                         ",
	"Numbers store a fixed point decimal number.                       ",
	"To insert a number, type in digits 0-9.                           ",
	"You can use a dot (.) or a comma (,) as the decimal separator;    ",
	"the decimal separator in the answer is a dot by default (.).      ",
	"Spaces are optional between numbers, signs and fractions. Numbers ",
	"support all operations listed in #Operators and many functions.   ",
	"                                                                  ",
	"# Vectors                                                         ",
	"A vector is an array that contains numbers. To enter a vector,    ",
	"type in numbers in parentheses, separated by semicolons (;).      ",
	"Vectors support addition, subtraction, scalar multiplication and a",
	"few named functions; for example ((2; 1) + (-1; 0.5)) = ((1; 1.5))",
	"                                                                  ",
	"# Unit values                                                     ",
	"In real life, we often deal with unit values rather than numbers. ",
	"To input a unit value, first enter its number part (optional), and",
	"then the unit name (220 V). The unit name can contain a SI prefix ",
	"(km, mA). All the SI units can be used, also the following        ",
	"non-decimal periods of time: (s), (min), (h), (day), (year). You  ",
	"can create your own units using (unit(\"name\"). Working with units ",
	"is similar to working with numbers, but some operations may cause ",
	"an error, for example (km + kg), (25 - J) or (11^day).            ",
	"                                                                  ",
	"# Angles                                                          ",
	"Theese objects contain angle values.                              ",
	"To type in an angle as an argument of a function, you can use     ",
	"radians or degrees; for example (sin 0.2618) and (sin 15°). The   ",
	"degree sign is entered with the '%' key. If you want to input an  ",
	"angle in radians outside of a function, use (angle(x)). Because   ",
	"angle is a special case of unit values, the same operations are   ",
	"available; angles can be arguments for trigonometric functions.   ",
	"                                                                  ",
	"# Strings                                                         ",
	"Strings store textual data. To type in a string, just write its   ",
	"content in double quotes (mass(\"H2O\")). Strings support addition, ",
	"but usually they are used as function inputs or outputs.          ",
	"                                                                  ",
	"# Operators                                                       ",
	"Operators perform simple operatoins, such as addition or division.",
	"The preference of the operators controls the order of evaluation: ",
	"operators with higher preference are calculated first. To specify ",
	"a particular order of evaluation, use parentheses. Exponentiation ",
	"is calculated right to left: (2^2^3) = (2^8). Sometimes, (*) sign ",
	"can be ommitted: (a√b), ((5+1)(5-1)) and (10 m) are equivalent to ",
	"(a * √b), ((5+1) * (5-1)) and (10 * m).                           ",
	"The following characters are entered by keystroke:                ",
	"{key}       {symbol}       {description}                          ",
	"+           (+)            plus sign                              ",
	"-           (-)            minus sign                             ",
	"*           (*)            multiplication sign                    ",
	":           (:)            division sign                          ",
	"^, **       (^),(**)       exponentiation sign                    ",
	"(, )        '(', ')'       parentheses                            ",
	"/                          horisontal fraction bar                ",
	"                                                                  ",
	"Preferences:                                                      ",
	"{operators}                {preference}                           ",
	"addition, subtraction      1                                      ",
	"multiplication, division   2                                      ",
	"named functions            3                                      ",
	"implicit multiplication    4                                      ",
	"unary minus                5                                      ",
	"exponentiation, roots      6                                      ",
	"fraction bar               --                                     ",
	"                                                                  ",
	"# Functions                                                       ",
	"Functions perform more complicated operations, and are typed in   ",
	"with the function name preceeding the arguments (fctr(5)). If     ",
	"there is more than one argument, they should be separated with    ",
	"semicolons (;); for example (sum(1; 2; 3)). Parentheses can be    ",
	"omitted, if it will not lead to an incorrect order of evaluation: ",
	"(ln 12 * 5) and (ln(12 * 5)) aren't equal; see Preferences in     ",
	"#Operators for more information. Trigonometric functions can be   ",
	"exponentiated: (tg^2 x) is equivalent to ((tg x)^2).              ",
	"The following functions can be entered as text:                   ",
	"{name}   {argument type}   {function}                             ",
	"sin         ang -> num     }                                      ",
	"cos         ang -> num     }                                      ",
	"tan, tg     ang -> num     } trigonometric and                    ",
	"cot, ctg    ang -> num     }                                      ",
	"arcsin      num -> ang     } inverse trigonometric                ",
	"arccos      num -> ang     }                                      ",
	"arctan/tg   num -> ang     }                                      ",
	"ln          num -> num     natural logorithm                      ",
	"abs         num -> num     absolute value                         ",
	"angle       num -> ang     convert number to angle in radians     ",
	"fctr        num -> num     factorial                              ",
	"pf          num -> str     prime factorisation                    ",
	"mass        str -> num     calculate chemical compound mass       ",
	"unit        str -> unit    create a unit                          ",
	"comb     num;num-> num     number of combinations, C(n; k)        ",
	"perm     num;num-> num     number of permutations, P(n; k)        ",
	"sum      vector -> num     sum of given numbers                   ",
	"average  vector -> num     arithmetic mean                        ",
	"norm     vector -> num     sum of squares of given numbers        ",
	"                                                                  ",
	"Note, that 'vector' means that the function accepts               ",
	"any number of arguments or a vector, containing them.             ",
	"                                                                  ",
	"# Commands                                                        ",
	"The following commands are entered by keystroke:                  ",
	"{key}       {command}                                             ",
	"?           help                                                  ",
	"Ctrl-x      exit                                                  ",
	"Ctrl-n      new equation                                          ",
	"Enter, =    show the answer                                       ",
	"Space       change answer notation    } only after Enter, =       ",
	"'           toggle greek alphabet letter input                    ",
	"                                                                  ",
	"Enter, /    go to the start of the numerator from denomenator     ",
	"Enter       go out of the current fraction                        ",
	"◂           go back one character                                 ",
	"▸           go forward one character                              ",
	"▴           go to the end of the numerator from denomenator       ",
	"▾           go to the end of the denomenator from numerator       ",
	"Ctrl-◂      go back one number                                    ",
	"Ctrl-▸      go forward one number                                 ",
	"Home        go to the beginning of the equation                   ",
	"End         go to the end of the equation                         ",
	"                                                                  ",
	"Bsp         delete the character to the left of the cursor        ",
	"Ctrl-Bsp    delete last number or fraction                        ",
	"Ctrl-p      parenthise all                                        ",
	"_           'forgotten' fraction bar                              ",
	"                                                                  ",
	"Ctrl-b      toggle debugging message (for developers)             ",
	"                                                                  ",
	"# Greek letter and rare symbol input                              ",
	"Sometimes you need to input values, represented by greek letters, ",
	"such as pi (π). To input a greek letter press ('), then input     ",
	"a letter of the english alphabet; for example ('+W) will give (Ω).",
	"Get the following symbols through (') + (given key):              ",
	"{key}: {symbol}                                                   ",
	"q: θ,    w: ω,    e: ε,    r: ρ,    t: τ,    m: μ,                ",
	"y: ψ,    u: υ,    i: ι,    o: ο,    p: π,    a: α,                ",
	"s: σ,    d: δ,    f: φ,    g: γ,    h: η,    k: κ,                ",
	"l: λ,    z: ζ,    x: χ,    c: ξ,    b: β,    n: ν,   ';': ς.      ",
	"                                                                  ",
	"Some other rare symbols and input shortcuts:                      ",
	"{key}       {symbols}      {description}                          ",
	"'+v         (√)            square root sign                       ",
	"'+V         (∛)            cube root sign                         ",
	"%           (°)            degree sign                            ",
	"Shift-2     (^2)           square                                 ",
	"Shift-3     (^3)           cube                                   ",
	"Shift-4     ( * 10^)       exponentiation                         ",
	"                                                                  ",
	"# Answer form                                                     ",
	"After pressing Enter or =, you will see the answer. The answer can",
	"be presented in several forms; to switch between them, use Space. ",
	"Here are the answer forms for numbers, vectors and units:         ",
	"{form}      {number / vector repr}    {unit value representation} ",
	"default     (x*10^3n), (x*10^n)       -practical-                 ",
	"percents    (x%)                      (x*10^n units)              ",
	"no-shift    (x)                       (x units)                   ",
	"                                                                  ",
	"The following are the answer forms for angles:                    ",
	"{form}      {symbol}       {answer representation}                ",
	"radians     (RAD)          angle in radians                       ",
	"degrees     (DEG)          angle in degrees, rounded to 3 digits  ",
	"pi/ang      (π:x)          a ratio between π radians and the angle",
	"                                                                  "]
	
	funcpad = [[' ?', 'help'], ['^X', 'exit'],
				['^N', 'new equation'], ['^P', 'parenthise all'],
				['Shift-2', 'x^2'], ['Shift-3', 'x^3'], 
				['Shift-4', 'x*10^y'], ['Shift-5', 'a°'],
				["'+v", '√x'], ["'", 'greek'],
				['/', 'fraction'], ['=', 'answer']]
	
	funcpad_greek = [['q:', 'θ'], ['w:', 'ω'], ['e:', 'ε'], ['r:', 'ρ'], ['t:', 'τ'],
			['y:', 'ψ'], ['u:', 'υ'], ['i:', 'ι'], ['o:', 'ο'], ['p:', 'π'], ['a:', 'α'],
			['s:', 'σ'], ['d:', 'δ'], ['f:', 'φ'], ['g:', 'γ'], ['h:', 'η'], ['k:', 'κ'],
			['l:', 'λ'], ['z:', 'ζ'], ['x:', 'χ'], ['c:', 'ξ'], ['b:', 'β'], ['n:', 'ν'],
			['m:', 'μ'], [';', 'ς'], ['', ''], ['v:', '√'], ['V:', '∛']]
	
	def __init__(self, stdscr):
		"""The initialiser of the class.
		
		Arguments:
		stdscr -- the terminal screen
		"""
		self.name = 'fractum'
		self.version = '1.4 '
		self.stln = ''
		
		# initialise the screen
		self.scr = stdscr
		self.line = 4
		self.col = 2
		self.new = True  # show welcome message
		self.debug = False # show debug message
		self.key = ''
		
		# initialise fractum's memory
		U = Unit.usr_init
		self.vars = {'x': Dec(10),      'Vm': Dec('22.4'),
			'day': U(Dec(86400), 's'),  'year': U(Dec(31557600), 's'),
			'c': Unit(Dec(299792458), {'m': 1, 's': -1})}
		self.vals = {'π': gl_pi, 'e': gl_e} | Unit.si_units()
		self.vals |= {'°': Angle.usr_init(1, 'degree')}
		
		# initialise the equation
		self.mode = 'normal'
		self.eqn = [' ']
		self.coordlist = [(0, self.col)]
		self.ansmode = 0
		
		# initialise the cursor
		curses.curs_set(2)
		self.cursor_num = 1
		self.cursor_tuple = (0, 1)
		
		# store rare symbols used in the calculator
		self.chars = AttributeSet()
		self.chars.bar = '─'     # opt: '-'
		self.chars.edge_l = '╶'  # opt: '─', '-', ' '
		self.chars.edge_r = '╴'  # opt: '─', '-', ' '
		
		# initialise keystroke input:
		self.chem_mode = False  # capital letters input
		# - keys to insert characters that aren't on the keyboard
		self.hotkeys = {'@': '^2', '#': '^3', '$': ' * 10^', '%': '°'}
		greek = 'ΑαΒβΓγΔδΕεΖζΗηΘθΙιΚκΛλΜμΝνΞξΟοΠπΡρΣσςΤτΥυΦφΧχΨψΩω'
		self.insert = {'q': 'θ', 'w': 'ω', 'e': 'ε', 'r': 'ρ', 't': 'τ',
			'y': 'ψ', 'u': 'υ', 'i': 'ι', 'o': 'ο', 'p': 'π', 'a': 'α', 
			's': 'σ', 'd': 'δ', 'f': 'φ', 'g': 'γ', 'h': 'η', 'k': 'κ',
			'l': 'λ', 'z': 'ζ', 'x': 'χ', 'c': 'ξ', 'b': 'β', 'n': 'ν',
			'm': 'μ'}
		capgreek = (lambda a: greek[greek.find(self.insert[a]) - 1])
		self.insert |= {t.upper(): capgreek(t) for t in self.insert}
		self.insert |= {'v': '√', 'V': '∛', ';': 'ς'}
		# - keys to insert directly
		self.keys = '0.,123456789()[]"*+-^:; '
		
		# initialise the calculator's theme
		self.txt_col = curses.COLOR_WHITE
		self.backgr_col = curses.COLOR_BLACK
		self.style_col = curses.COLOR_GREEN
		self.error_col = curses.COLOR_RED
		# - title and statusline color
		curses.init_pair(1, self.backgr_col, self.style_col)
		# - equation color
		curses.init_pair(2, self.txt_col, self.backgr_col)
		# - answer color
		curses.init_pair(3, self.style_col, self.backgr_col)
		# - error color
		curses.init_pair(4, self.error_col, self.backgr_col)
		
		# initialise the help text
		self.helpline = 0
		self.helppad = curses.newpad(len(Calculator.helptext), 70)
		self.write_help()
	
	def quit(self):
		"""Exit the calculator."""
		self.scr.clear()
		self.scr.refresh()
		sys.exit()
	
	def write_help(self):
		self.helppad.attron(curses.color_pair(2))
		def get_col_pair(string):
			if string.startswith('# '):
				return curses.color_pair(3)
			if string.startswith('{'):
				return curses.color_pair(3)
			return curses.color_pair(2)
		line = 0
		for string in Calculator.helptext:
			self.helppad.addstr(line, 1, string, get_col_pair(string))
			line += 1
	
	def write_statusline(self):
		"""Write the statusline string according to the situation."""
		self.stln = ''
		if self.mode == 'normal':
			if self.new:
				self.stln = 'Welcome to fractum! Press (?) for help'
			elif self.eqn == [' ']:
				self.stln = 'Write an equation'
			elif self.cursor_tuple[0] % 3 == 1:  # in numerator
				self.stln = 'Enter, Down: go to denomenator'
			elif self.cursor_tuple[0] % 3 == 2:  # in denomenator
				self.stln = 'Enter: go out of fraction,  Up: go to numerator'
			else:
				self.stln = 'Enter: see the answer'
		elif self.mode == 'insert':
			self.stln = 'Key: greek letter,  Shift-key: capital greek letter  (hints below)'
		elif self.mode == 'help':
			self.stln = 'Up/Down: scroll,  any key: get back'
		elif self.mode == 'final':
			self.stln = 'Backspace: get back,  Space: answer form,  Enter: new,  any key: exit'
	
	def usr_input(self):
		"""Wait for the user to press a key and analise it."""
		key = self.scr.getkey()
		if self.command_read(key):
			pass
		elif self.mode == 'normal':
			self.arrows_read(key)
			self.normal_mode_read(key)
		elif self.mode == 'insert':
			self.insert_mode_read(key)
		elif self.mode == 'help':
			self.help_mode_read(key)
		elif self.mode == 'final':
			self.final_mode_read(key)
	
	def eqn_addstr(self, a, b, string):
		"""Insert a string in the calculator's equation.
		
		Arguments:
		a, b -- position of insertion
		string -- the string to insert
		"""
		if len(self.eqn[a]) == b-1:
			self.eqn[a] += string
		else:
			self.eqn[a] = self.eqn[a][:b-1] + string + self.eqn[a][b-1:]
		self.cursor_num += len(string)
	
	def command_read(self, key):
		"""Run command by key."""
		if len(key) != 1:
			return None
		if key == curses.ascii.ctrl('X'):
			self.quit()
		elif key == '?':
			self.mode = 'help'
			curses.curs_set(0)
		elif key == curses.ascii.ctrl('P'):
			a, b = self.cursor_tuple
			self.eqn[0] = '(' + self.eqn[0]
			self.eqn[-1] = self.eqn[-1][:-1] + ') '
			if a == 0 and b == 1:
				pass
			elif a == len(self.eqn) - 1 and b == len(self.eqn[a]) - 2:
				self.cursor_num += 2
			else:
				self.cursor_num += 1
		elif key == curses.ascii.ctrl('B'):
			self.debug = not self.debug
		elif key == curses.ascii.ctrl('N'):
			ans = self.ans_calc()
			num_types = (int, Angle, Dec, Unit)
			if self.mode == 'final' and type(ans) in num_types:
				self.vars |= {'ans': ans}
			self.mode = 'normal'
			self.ansmode = 0
			self.eqn = [' ']
			self.coordlist = [(0, self.col)]
		else:
			return False
		return True
	
	def normal_mode_read(self, key):
		"""Change the equation according to key for normal mode."""
		a, b = self.cursor_tuple
		if key in self.keys or key.isalpha():
			self.eqn_addstr(a, b, key)
		elif key in list(self.hotkeys):
			self.eqn_addstr(a, b, self.hotkeys[key])
		elif key == '/':
			if a % 3 == 0:
				eqn2 = self.eqn[:a] + [self.eqn[a][:b-1] + ' ']
				eqn2 += [' ', ' '] + [self.eqn[a][b-1:]] + self.eqn[a+1:]
				self.eqn = eqn2
				self.cursor_num += 1
			if a % 3 == 1:
				self.cursor_num += len(self.eqn[a]) - b + 1
		elif key == '_':
			self.underscore()
		elif key == "'":
			self.mode = 'insert'
		elif key in '\n\\':
			if a == len(self.eqn) - 1:
				self.mode = 'final'
			elif self.eqn[a + 1] == ' ':
				self.cursor_num += len(self.eqn[a]) - b + 1
			elif self.eqn[a + 2] == ' ':
				self.cursor_num += len(self.eqn[a]) + len(self.eqn[a+1]) - b + 1
			else:
				self.mode = 'final'
		elif key == 'KEY_BACKSPACE':
			self.backspace(a, b)
		elif key == chr(8):  # Ctrl + KEY_BACKSPACE
			self.ctrl_backspace(a, b)
		elif key == '=':
			self.mode = 'final'
		self.key = key
	
	def arrows_read(self, key):
		"""Change cursor position according to key for normal mode."""
		a, b = self.cursor_tuple
		if key == 'KEY_RIGHT':
			self.cursor_num += 1
		elif key == 'KEY_LEFT':
			self.cursor_num -= 1
		elif key == 'KEY_UP':
			if a % 3 == 2:
				self.cursor_num -= b
		elif key == 'KEY_DOWN':
			if a % 3 == 1:
				self.cursor_num += len(self.eqn[a])
				self.cursor_num += len(self.eqn[a+1]) - b
		elif key == 'KEY_HOME':
			self.cursor_num = 1
		elif key == 'KEY_END':
			self.cursor_num = sum([len(s) for s in self.eqn])
		elif key == 'kRIT5':  # Ctrl + KEY_RIGHT
			self.ctrl_right_arrow(a, b)
		elif key == 'kLFT5':  # Ctrl + KEY_LEFT
			self.ctrl_left_arrow(a, b)
	
	def ctrl_right_arrow(self, a, b):
		"""Change cursor position a number to the right.
		
		Arguments:
		a, b -- initial position of cursor
		"""
		def skipword(smbs):
			if self.eqn[a][b-1] in smbs:
				p = 0
				while self.eqn[a][b+p-1] in smbs:
					p += 1
				self.cursor_num += p - 1
		
		if a % 3 == 0:
			if b == len(self.eqn[a]):
				if a != len(self.eqn) - 1:
					self.cursor_num += len(self.eqn[a+1])
					self.cursor_num += len(self.eqn[a+2])
			elif self.eqn[a][b-1] == ' ':
				p = 0
				while self.eqn[a][b+p-1] == ' ':
					p += 1
				self.cursor_num += p
				b += p
			skipword('1234567890,.')
			skipword('+-*:()^[]""')
			skipword('abcdefghijklmnopqrstuvwxyz')
			skipword('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
		elif a % 3 == 1:
			self.cursor_num += len(self.eqn[a])
			self.cursor_num += len(self.eqn[a+1]) + b
		else:
			self.cursor_num += len(self.eqn[a]) - b
		self.cursor_num += 1
	
	def ctrl_left_arrow(self, a, b):
		"""Change cursor position a number to the left.
		
		Arguments:
		a, b -- initial position of cursor
		"""
		def skipword(smbs):
			nonlocal b
			if self.eqn[a][b-2] in smbs:
				p = 0
				while p < b and self.eqn[a][b-p-2] in smbs:
					p += 1
				self.cursor_num -= p - 1
				b -= p
		
		if a % 3 == 0:
			if b == 1 and a != 0:
				self.cursor_num -= len(self.eqn[a-1])
				self.cursor_num -= len(self.eqn[a-2])
			skipword('1234567890,.')
			skipword('+-*:()^[]""')
			skipword('abcdefghijklmnopqrstuvwxyz')
			skipword('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
			if self.eqn[a][b-2] == ' ':
				p = 0
				while p < b and self.eqn[a][b-p-2] == ' ':
					p += 1
				self.cursor_num -= p
		elif a % 3 == 1:
			self.cursor_num -= b - 1
		else:
			self.cursor_num -= len(self.eqn[a-1]) + b - 1
		self.cursor_num -= 1
	
	def underscore(self):
		"""Insert a 'forgotten' fraction bar in the equation."""
		if len(self.eqn) % 3 == 1:
			for i in range(len(self.eqn[-1])):
				if self.eqn[-1][i] not in ' ^*:+-':
					part1 = self.eqn[-1][:i] + ' '
					part2 = self.eqn[-1][i:]
					self.eqn = self.eqn[:-1] + [part1, part2, ' ', ' ']
					self.cursor_num += 2
					break
	
	def backspace(self, a, b):
		"""Delete the character before cursor."""
		if b == 1:
			# don't delete anything if cursor is on the 1-st char
			if a == 0:
				return None
			# delete an empty fraction if cursor is in front of it
			elif a % 3 == 0:
				if self.eqn[a] == ' ' and self.eqn[a-1] == ' ':
					center = self.eqn[a-3][:-1] + self.eqn[a][:-1] + ' '
					self.eqn = self.eqn[:a-3] + [center] + self.eqn[a+1:]
			# delete an empty fraction from it's numerator
			elif a % 3 == 1:
				if self.eqn[a+1] == ' ':
					center = self.eqn[a-1][:-1] + self.eqn[a][:-1]
					center += self.eqn[a+2]
					self.eqn = self.eqn[:a-1] + [center] + self.eqn[a+3:]
			self.cursor_num -= 1
			return None
		else:
			self.eqn[a] = self.eqn[a][:b-2] + self.eqn[a][b-1:]
			self.cursor_num -= 1
	
	def ctrl_backspace(self, a, b):
		"""Delete the number before cursor."""
		
		def delword(smbs):
			if self.eqn[a][b-2] in smbs:
				p = 0
				while p < b and self.eqn[a][b-p-2] in smbs:
					p += 1
				while p < b-1 and self.eqn[a][b-p-2] == ' ':
					p += 1
				if b == len(self.eqn[a]):
					part = [self.eqn[a][:b-p-1] + self.eqn[a][b:] + ' ']
					self.eqn = self.eqn[:a] + part + self.eqn[a+1:]
					self.cursor_num -= p - 1
				else:
					part = [self.eqn[a][:b-p-1] + ' ' + self.eqn[a][b:]]
					self.eqn = self.eqn[:a] + part + self.eqn[a+1:]
					self.cursor_num -= p
		
		if b == 1:
			# don't delete anything if cursor is on the 1-st char
			if a == 0:
				return None
			# delete a fraction if cursor is in front of it
			elif a % 3 == 0:
				self.cursor_num -= len(self.eqn[a-2])
				self.cursor_num -= len(self.eqn[a-1])
				center = self.eqn[a-3][:-1] + self.eqn[a]
				self.eqn = self.eqn[:a-3] + [center] + self.eqn[a+1:]
			# delete an empty fraction from it's numerator
			elif a % 3 == 1:
				if self.eqn[a+1] == ' ':
					center = self.eqn[a-1][:-1] + self.eqn[a][:-1]
					center += self.eqn[a+2]
					self.eqn = self.eqn[:a-1] + [center] + self.eqn[a+3:]
			# delete a fraction bar from the denomenator
			elif a % 3 == 2:
				center = self.eqn[a-2][:-1] + self.eqn[a-1][:-1]
				center += self.eqn[a]
				self.eqn = self.eqn[:a-2] + [center] + self.eqn[a+2:]
			self.cursor_num -= 1
			return None
		else:
			if self.eqn[a][b-2] in '0123456789.,':
				delword('0123456789.,')
			elif self.eqn[a][b-2].islower():
				delword('abcdefghijklmnopqrstuvwxyz')
			elif self.eqn[a][b-2].isupper():
				delword('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
			else:
				self.backspace(a, b)
				b -= 1
				if self.eqn[a][b-2] == ' ':
					p = 0
					while p < b and self.eqn[a][b-p-2] == ' ':
						p += 1
					self.cursor_num -= p
			if a % 3 != 0:
				self.cursor_num -= 1
	
	def insert_mode_read(self, key):
		"""Insert a symbol according to key."""
		a, b = self.cursor_tuple
		if key in self.insert:
			self.eqn_addstr(a, b, self.insert[key])
		elif key == 'KEY_RESIZE':
			return None
		self.mode = 'normal'
	
	def help_mode_read(self, key):
		"""Change position of the help text according to key."""
		maxy, maxx = self.scr.getmaxyx()
		if key == 'KEY_UP':
			if self.helpline > 0:
				self.helpline -= 1
		elif key == 'KEY_DOWN':
			if self.helpline < len(Calculator.helptext) - maxy + 5:
				self.helpline += 1
		elif key != 'KEY_RESIZE':
			self.mode = 'normal'
			curses.curs_set(2)
	
	def final_mode_read(self, key):
		"""Exit or return to normal mode according to key."""
		if key in ['KEY_BACKSPACE', 'KEY_LEFT']:
			self.mode = 'normal'
			self.ansmode = 0
		elif key == ' ':
			self.ansmode = (self.ansmode + 1) % 3
		elif key == '\n':
			self.command_read(curses.ascii.ctrl('N'))
		elif key != 'KEY_RESIZE':
			self.quit()
	
	def cursor_update(self):
		"""Change the list cursor position according to it's number position."""
		eq_len = sum([len(s) for s in self.eqn])
		if self.cursor_num < 1:
			self.cursor_num = 1
		elif self.cursor_num > eq_len:
			self.cursor_num = eq_len
		c = 0
		a = -1
		for eq_str in self.eqn:
			a += 1
			if c + len(eq_str) >= self.cursor_num:
				self.cursor_tuple = (a, self.cursor_num - c)
				return None
			c += len(eq_str)
	
	def del_last_space(self, data):
		"""Remove the last characters of strings in data.
		
		If data is a string, return it with the last character removed.
		"""
		if type(data) is str:
			return data[:-1]
		else:
			return [string[:-1] for string in data]
	
	def printeq(self, y, x, msg, color=2):
		"""Print a string on the screen.
		
		Arguments:
		y, x -- position according to the equation start position
		msg -- a string to print
		color -- number of text's color_pair (default=2)
		"""
		self.scr.addstr(self.line + y, x, msg, curses.color_pair(color))
	
	def coordlist_update(self):
		"""Change the positions of the equation pieces."""
		p = self.del_last_space
		coordlist = []
		x = self.col
		for b in range(len(self.eqn) // 3):
			a = b * 3
			mainstr, numstr, denomstr = p(self.eqn[a:a+3])
			coordlist.append((0, x))
			x += len(mainstr)
			fractlen = max(len(numstr), len(denomstr)) + 2
			numx = x + math.ceil((fractlen - len(numstr) - 1) / 2)
			denomx = x + math.ceil((fractlen - len(denomstr) - 1) / 2)
			coordlist.append((-1, numx))
			coordlist.append((1, denomx))
			x += fractlen
		if len(self.eqn) % 3 == 1:
			coordlist.append((0, x))
		elif len(self.eqn) % 3 == 2:
			mainstr, numstr = p(self.eqn[-2:])
			coordlist.append((0, x))
			x += len(mainstr)
			coordlist.append((-1, x + 1))
		self.coordlist = coordlist
	
	def scr_update(self):
		"""Update the calculator's screen."""
		self.scr.clear()
		self.print_title()
		self.print_statusline()
		self.print_functions()
		if self.mode == 'help':
			self.scr.refresh()
			self.print_help()
			return None
		self.print_answer_mode()
		if self.debug:
			self.print_debug()  # HELP
		try:
			self.print_equation()
			if self.mode == 'final':
				self.print_answer()
				self.scr.refresh()
				return None
			self.move_cursor()
		except Exception:
			self.print_error()
			raise  # DEBUG
			return None
		self.scr.refresh()
	
	def print_title(self):
		"""Print the calculator's title line on the screen."""
		maxy, maxx = self.scr.getmaxyx()
		name = self.name + ' ' + self.version
		space = (maxx - len(name)) // 2
		last = ' '
		if maxx % 2 == len(name) % 2:
			last = ''
		title = ' ' * space + name + ' ' * space + last
		self.scr.addstr(0, 0, title, curses.color_pair(1))
	
	def print_statusline(self):
		"""Print the statusline on the screen."""
		maxy, maxx = self.scr.getmaxyx()
		self.write_statusline()
		space = maxx - len(self.stln) - 1
		stln = ' ' + self.stln + ' ' * space
		self.scr.addstr(maxy - 3, 0, stln, curses.color_pair(1))
	
	def print_functions(self):
		"""Print the keyboard shortcuts on the screen."""
		maxy, maxx = self.scr.getmaxyx()
		rows = 2
		funcpad = self.funcpad
		if self.mode == 'insert':
			funcpad = self.funcpad_greek
		cols = (len(funcpad) + (rows - 1)) // rows + 1
		total = 0
		for x in range(cols - 1):
			lens = []
			for y in range(rows):
				leny = len(funcpad[rows * x + y][0])
				leny += len(funcpad[rows * x + y][1])
				lens.append(leny)
			total += max(lens)
		space = (maxx - total) // cols
		if space < 3:
			return None
		x = (maxx - space * cols - total) // 2 + space
		for s in range(cols - rows + 1):
			lens = []
			for t in range(rows):
				lens.append(self.print_func(maxy - rows + t, x, s * rows + t, funcpad))
			x += max(lens) + space
		try:
			s += 1
			for t in range(rows):
				lens.append(self.print_func(maxy - rows + t, x, s * rows + t, funcpad))
		except IndexError:
			pass
	
	def print_func(self, y, x, n, funcpad):
		"""Print one keyboard shortcut on the screen.
		
		Arguments:
		y, x -- position of the shortcut on the screen
		n -- number of shortcut
		shift -- start of the shortcut
		"""
		pf = funcpad[n]
		self.scr.addstr(y, x - 1, pf[0], curses.color_pair(3))
		self.scr.addstr(y, x + len(pf[0]), pf[1], curses.color_pair(2))
		return len(pf[0]) + len(pf[1])
	
	def print_answer_mode(self):
		"""Print the answer mode switcher on the screen."""
		if self.mode != 'final':
			self.scr.addstr(1, 5, 'equation', curses.color_pair(3))
			return None
		def col_pair(n):
			if self.ansmode == n:
				return curses.color_pair(1)
			return curses.color_pair(3)
		ans = self.ans_calc()
		if type(ans) == str:
			if ans.startswith('$'):
				self.scr.addstr(1, 5, ' ERROR ', curses.color_pair(1))
				return None
			anstype = 'string'
			modes = ['', '', '']
		elif type(ans) == np.ndarray:
			anstype = 'vector'
			modes = [' DECIMAL ', ' PERCENTS ', ' SIMPLE ']
		elif type(ans) == Dec:
			anstype = 'number'
			modes = [' DECIMAL ', ' PERCENTS ', ' SIMPLE ']
		elif type(ans) == Unit:
			anstype = 'unit value'
			modes = [' PREFIX ', 'NO-PREFIX', ' SIMPLE ']
		elif type(ans) == Angle:
			anstype = 'angle'
			modes = [' DEG ', ' RAD ', ' π : RAD ']
		else:
			return None
		
		self.scr.addstr(1, 5, anstype, curses.color_pair(3))
		x = 20
		for a in range(3):
			self.scr.addstr(1, x, modes[a], col_pair(a))
			x += len(modes[a]) + 5
	
	def print_debug(self):
		"""Print some extra info about the calculator's current state."""
		maxy, maxx = self.scr.getmaxyx()
		try:
			self.scr.addstr(maxy - 7, 2, str(self.simple_notation()))
			self.scr.addstr(maxy - 6, 2, str(self.infix_notation()))
			self.scr.addstr(maxy - 5, 2, str(self.postfix_notation()))
			k = self.key
			if self.key == '\n':
				k = 'RETURN'
			self.scr.addstr(maxy - 8, 2, str(f'{k} {ord(self.key)}'))
		except Exception:
			pass
	
	def print_equation(self):
		"""Print the equation on the calculator screen."""
		p = self.del_last_space
		for s in range(len(self.eqn)):
			y, x = self.coordlist[s]
			self.printeq(y, x, p(self.eqn[s]))
			if y == 0 and s < len(self.coordlist) - 2:
				self.print_fract_line(s)
	
	def print_fract_line(self, s):
		"""Print the horisontal fraction bar with index s."""
		p = self.del_last_space
		y, x = self.coordlist[s]
		mainstr, numstr, denomstr = p(self.eqn[s:s+3])
		fractlen = max(len(numstr), len(denomstr)) + 2
		middle = self.chars.bar * (fractlen - 2)
		bar = self.chars.edge_l + middle + self.chars.edge_r
		if numstr == '' and denomstr == '':
			bar = self.chars.edge_l + self.chars.bar
		self.printeq(y, x + len(mainstr), bar)
	
	def print_help(self):
		"""Show the help text."""
		maxy, maxx = self.scr.getmaxyx()
		self.helppad.refresh(self.helpline, 0, 1, 0, maxy - 5, 70)
	
	def print_answer(self):
		"""Print the answer on the screen."""
		p = self.del_last_space
		y, x = self.coordlist[-1]
		x += len(p(self.eqn[-1]))
		ans = self.ans_format(self.ans_calc())
		self.printeq(0, x, ' = ')
		x += 3
		if type(ans) is tuple:
			self.printeq(-1, x, ans[0], color=3)
			self.printeq(1, x, ans[2], color=3)
			self.printeq(0, x, ans[1], color=3)
		elif ans.startswith('$'):
			self.printeq(0, x, ans[1:], color=4)
		else:
			self.printeq(0, x, ans, color=3)
	
	def print_error(self):
		"""Print the display error text."""
		self.scr.clear()
		msg = 'DISPLAY ERROR. Try one of theese actions:\n '
		msg += 'Bsp: delete last character, ^X: exit,  ^N: restart.\n'
		msg += 'Making the screen bigger may also help, if possible.'
		self.scr.addstr(msg, curses.color_pair(4))
		self.scr.refresh()
	
	def move_cursor(self):
		"""Move the cursor to it's assigned position."""
		a, b = self.cursor_tuple
		x = self.coordlist[a][1] + b - 1
		y = self.coordlist[a][0] + self.line
		curses.curs_set(2)
		if self.mode == 'insert':
			#curses.curs_set(0)
			self.scr.addstr("'")
		self.scr.move(y, x)
	
	def check_names(self):
		"""Check the equation for unknown names."""
		words = re.findall('[a-z]+', ''.join(self.eqn))
		for word in words:
			if word not in list(Calculator.operdict | self.vars | self.vals):
				return word
		return ''
	
	def check_parentheses(self):
		"""Check the equation for unmatched parenthesses."""
		ps = re.findall('[()]', ''.join(self.eqn))
		stack = 0
		for p in ps:
			if p == ')':
				if stack == 0:
					return False
				else:
					stack -= 1
			elif p == '(':
				stack += 1
		return stack == 0
	
	def check_sq_brackets(self):
		"""Check the equation for unmatched square brackets."""
		ps = re.findall(r'[\[\]]', ''.join(self.eqn))
		stack = 0
		for p in ps:
			if p == ']':
				if stack == 0:
					return False
				else:
					stack -= 1
			elif p == '[':
				stack += 1
		return stack == 0
	
	def simple_notation(self):
		"""Transform the equation into a string.
		
		Output syntax:
		+ - * : -- arithmetic operators
		. -- decimal separator
		^ -- exponentiation sign
		@ -- default multiplication sign
		~ -- unary minus
		√ ∛ -- roots
		() -- parentheses
		'string' -- string
		name -- name
		"""
		line = ''
		for b in range(len(self.eqn) // 3):
			a = b * 3
			line += self.eqn[a][:-1] + '( '
			line += self.eqn[a+1][:-1] + ' )/( ' + self.eqn[a+2][:-1] + ' )'
		line += self.eqn[-1]
		# unificate all argument separators
		line = re.sub(r', ', r' ; ', line)
		line = re.sub(r';', r' ; ', line)
		# unificate all decimal separators
		line = re.sub(r',', r'.', line)
		# unificate all exponentiation symbols
		line = re.sub(r'\*\*', r'^', line)
		# delete unary '+' at start of the string
		line = re.sub(r'\A\s*\+', r'', line)
		# delete unary plus after an operator symbol
		line = re.sub(r'(['+Calculator.oper_regexp+r'(;])\s*\+', r'\1', line)
		# replace unary '-' at start of the string
		line = re.sub(r'\A\s*-', r'~', line)
		# replace unary '-' after an operator symbol
		line = re.sub(r'(['+Calculator.oper_regexp+r'(;])\s*-', r'\1~', line)
		# find all strings
		strings = re.findall(r'"[^"]*"', line)
		line = re.sub(r'"[^"]*"', '%', line)
		# find all function exponentiation
		for f in ['sin', 'cos', 'tg', 'ctg', 'tan', 'cot']:
			line = re.sub(f+r'\s*\^\s*([0-9.]+)', f+r'_ \1 |', line)
			for v in list(self.vars | self.vals):
				if v not in line:
					continue
				line = re.sub(f+r'\s*\^\s*('+v+')', f+r'_ \1 |', line)
			line = re.sub(f+r'\s*\^', f+'_ ', line)
		# work with names:
		for v in list(self.vars | self.vals):
			if v not in line:
				continue
			orx = '([' + Calculator.oper_regexp + r'\(\)])'
			# insert default multiplication symbol (@) between numbers and names
			line = re.sub(r'([0-9.]+)\s*'+v+r'([^a-zA-Z])', r'\1@'+v+r'\2', line)
			# insert @ between a name and a root symbol
			line = re.sub(r'([^a-zA-Z])'+v+r'\s*([√∛])', r'\1'+v+r'@\2', line)
			# insert @ between names and brackets
			line = re.sub(r'([^a-zA-Z])'+v+r'\s*\(', r'\1'+v+r'@(', line)
			line = re.sub(r'\)\s*'+v+r'([^a-zA-Z])', r')@'+v+r'\1', line)
			# insert @ between two names
			for v2 in list(self.vars | self.vals):
				if v2 not in line:
					continue
				line = re.sub('([^a-zA-Z])'+v+r'\s+'+v2+'([^a-zA-Z])', r'\1'+v+'@'+v2+r'\2', line)
		# insert @ between a number and a root symbol
		line = re.sub(r'([0-9.]+)\s*([√∛])', r'\1@\2', line)
		# insert @ between numbers and brackets
		line = re.sub(r'([0-9.]+)\s*\(', r'\1@(', line)
		line = re.sub(r'\)\s*([0-9.]+)', r')@\1', line)
		# insert @ between brackets
		line = re.sub(r'\)\s*\(', r')@(', line)
		# omit all '|'
		line = re.sub(r'\|', r' ', line)
		# unificate the number of spaces
		line = re.sub(r'(['+Calculator.oper_regexp+r'()])', r' \1 ', line)
		line = re.sub(r'\s+', r' ', line)
		# insert omitted strings
		for s in strings:
			line = re.sub(r'%', " '" + s[1:-1] + "'", line, count=1)
		return line
	
	def infix_notation(self):
		"""Transform the simple notation string into a list of tokens."""
		line = self.simple_notation()
		ls = re.split(r' ', line)
		d = 0
		for n in range(len(ls)):
			if ls[n-d] in ['', ' ']:
				del ls[n-d]
				d += 1
			elif re.fullmatch(r'[0-9.]+', ls[n-d]):
				if re.fullmatch(r'\.', ls[n-d]) or len(re.findall(r'\.', ls[n-d])) > 1:
					raise ValueError('incorrect number notation')
				ls[n-d] = Dec(ls[n-d])
			elif ls[n-d][0] + ls[n-d][-1] == "''":
				ls[n-d] = ls[n-d][1:-1]
			elif ls[n-d] in list(Calculator.operdict):
				ls[n-d] = Calculator.operdict[ls[n-d]]
			elif ls[n-d] in list(self.vals):
				ls[n-d] = self.vals[ls[n-d]]
			elif ls[n-d] in list(self.vars):
				ls[n-d] = self.vars[ls[n-d]]
			elif ls[n-d] not in '()':
				raise ValueError('unknown name: ' + ls[n-d])
		return ls
	
	def postfix_notation(self):
		"""Transform the infix notation into postfix notation."""
		ls = self.infix_notation()
		oper_stack = []
		output = []
		argsbreak = False
		for token in ls:
			if type(token) in [Dec, Angle, Unit]:
				output.append(token)
			elif type(token) is Operator:
				if token.left:
					while oper_stack and oper_stack[-1] != '(' and \
				token.pref < oper_stack[-1].pref:
						output.append(oper_stack.pop())
					oper_stack.append(token)
				else:
					while oper_stack and oper_stack[-1] != '(' and \
				token.pref <= oper_stack[-1].pref:
						output.append(oper_stack.pop())
					oper_stack.append(token)
			elif token == '(':
				oper_stack.append(token)
			elif token == ')':
				while oper_stack[-1] != '(':
					output.append(oper_stack.pop())
				oper_stack.pop()
			elif type(token) is str:
				output.append(token)
		return output + oper_stack[::-1]
	
	def ans_calc(self):
		"""Calculate the answer of current equation."""
		line = self.simple_notation()
		if not self.check_parentheses():
			return '$unmatched parentheses'
		if not self.check_sq_brackets():
			return '$unmatched square brackets'
		try:
			ls = self.postfix_notation()
			data_stack = []
			semicolons = 0
			for token in ls:
				if token == Operator.SEMICOLON:
					semicolons += 1
				else:
					if semicolons != 0:
						if len(data_stack) < semicolons + 1:
							return '$compilation error'
						g = []
						for _ in range(semicolons + 1):
							a = data_stack.pop()
							g.insert(0, a)
						data_stack.append(np.array(g))
						semicolons = 0
					if type(token) in [Dec, Angle, Unit, str]:
						data_stack.append(token)
					elif type(token) is Operator:
						if token.arg_num == 1:
							if len(data_stack) < 1:
								return '$compilation error'
							x = data_stack.pop()
							data_stack.append(token.calc(x))
						elif token.arg_num == 2:
							if len(data_stack) < 2:
								return '$compilation error'
							b = data_stack.pop()
							a = data_stack.pop()
							data_stack.append(token.calc(a, b))
			if len(data_stack) == 1:
				return data_stack[0]
			elif len(data_stack) == 0:
				return '$none'
			else:
				return np.array(data_stack)
		except SyntaxError as ex:
			return '$' + str(ex)
		except TypeError as ex:
			return '$' + str(ex)
		except ValueError as ex:
			return '$' + str(ex)
		except Unit.InvalidOperationError as ex:
			return '$' + str(ex)
		except ZeroDivisionError as ex:
			return '$division by zero: ' + str(ex)
		except decimal.InvalidOperation:
			return '$function domain error'
		except decimal.Overflow:
			return '$large number'
	
	def ans_format(self, ans):
		"""Transform ans to the required format."""
		if type(ans) is Angle:
			if ans.getpow('ang') != 1:
				part = '$incorrect answer type: angle^'
				ans = part + str(ans.getpow('ang'))
			elif self.ansmode == 0:  # 'DEG'
				ans = str(ans.degree()) + '°'
			elif self.ansmode == 1:  # 'RAD'
				ans = str(ans)
			elif self.ansmode == 2:  # 'PIRAD'
				ans = 'π : ' + str(ans.pirad())
		elif type(ans) is Unit:
			return ans.smart_str(self.ansmode)
		elif type(ans) is decimal.Decimal:
			if self.ansmode == 0:
				ans = numberform(ans)
			elif self.ansmode == 1:
				ans = str(round(ans * 100, gl_prec)) + '%'
			elif self.ansmode == 2:
				ans = str(remove_exponent(ans))
		elif type(ans) is np.ndarray:
			if self.ansmode == 0:
				ansformat = numberform
			elif self.ansmode == 1:
				ansformat = (lambda a: str(round(a * 100, gl_prec)) + '%')
			elif self.ansmode == 2:
				ansformat = (lambda a: str(remove_exponent(a)))
			g = ''
			for x in ans:
				if type(x) is Dec:
					g += ansformat(x) + '; '
				elif type(x) is Unit:
					g += x.smart_str() + '; '
				else:
					return '$incorrect answer type in array: ' + str(type(ans))
			return '(' + g[:-2] + ')'
		elif type(ans) is str:
			return ans
		else:
			ans = '$incorrect answer type: ' + str(type(ans))
		return ans
	
	def mainloop(self):
		"""Interact with the user."""
		while True:
			self.coordlist_update()
			self.cursor_update()
			self.scr_update()
			self.usr_input()
			self.new = False


wrapper(lambda stdscr: Calculator(stdscr).mainloop())
