import sys

ALIGN_RIGHT = 'r'
ALIGN_LEFT = 'l'
ALIGN_CENTER = 'c'

class TablePrinter( object ):
	_ROW_PREFIX = '|'
	_COLUMN_PREFIX = ' '
	_COLUMN_SUFFIX = ' |'
	_CORNER = '+'
	_LINE_CHAR = '-'

	def __init__( self, rows, headRow = None, alignmentsString = None, stream = sys.stdout ):
		if headRow is not None:
			headRow = [ str( thing ) for thing in headRow ]
		self._head = headRow
		self._rows = self._makeStringsFromRows( rows ) 
		self._columnWidths = []
		self._stream = stream
		if alignmentsString is None:
			alignmentsString = ALIGN_LEFT * len( self._rows[ 0 ] )
		self._alignments = list( alignmentsString )

	def _makeStringsFromRows( self, rows ):
		newRows = []
		for row in rows:
			newRows.append( [ str( x ) for x in row ] )
		return newRows

	def _findColumnWidths( self ):
		rows = self._rows
		if self._head:
			rows = [ self._head ] + rows
		self._columnWidths = [ 0 ] * len( rows[ 0 ] )
		for row in rows:
			widths = [ len( str( word ) ) for word in row ]
			self._columnWidths = [ max( self._columnWidths[ i ], widths[ i ] ) for i in range( len( widths ) ) ]
		self._width = self._calculateWidth()

	def output( self, flush = False ):
		self._outputTable()
		if flush:
			self._stream.flush()

	def _outputTable( self ):
		self._findColumnWidths()
		self._outputHead()
		for row in self._rows:
			self._outputRow( row )
		self._outputTail()

	def _outputRow( self, row ):
		line = self._ROW_PREFIX
		for i in range( len( row ) ):
			word = self._align( row[ i ], self._columnWidths[ i ], self._alignments[ i ] )
			line += self._COLUMN_PREFIX + word + self._COLUMN_SUFFIX
		self._writeLine( line )

	def _align( self, word, width, alignment ):
		words = { 
			ALIGN_LEFT:		word.ljust( width ),
			ALIGN_RIGHT: 	word.rjust( width ),
			ALIGN_CENTER: 	word.center( width ) }
		return words[ alignment ]

	def _outputHead( self ):
		self._horizontalLine( corners = True )
		if self._head:
			self._outputRow( self._head )
			self._horizontalLine( corners = True )

	def _outputTail( self ):
		self._horizontalLine( corners = True )

	def _calculateWidth( self ):
		return len( self._ROW_PREFIX ) + \
						 len( self._columnWidths ) * len( self._COLUMN_PREFIX ) + \
						 sum( self._columnWidths ) + \
						 len( self._columnWidths ) * len( self._COLUMN_SUFFIX )

	def _writeLine( self, line ):
		self._stream.write( '%s\n' % line )

	def _horizontalLine( self, corners = False ):
		if corners:
			line = self._CORNER + self._LINE_CHAR * ( self._width - 2 ) + self._CORNER
		else:
			line = self._LINE_CHAR * self._width
		self._writeLine( line )

if __name__ == '__main__':
	head = [ 'sex', 'name', 'height' ]
	rows = [ 
			[ 'M', 'Moshe Zuchmir', 167 ],
			[ 'F', 'Zoot Dingo',	208 ],
			[ 'F', 'Vera Xenon', 25 ]
										]

	print "with head:"
	printer = TablePrinter( rows, head, 'llr' )
	printer.output()

	print ''
	print "without head:"
	printer = TablePrinter( rows, alignmentsString = 'llr' )
	printer.output()
