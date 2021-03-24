//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

package explicit.result;

public class ProcessingResult<T>
{
	public static final double PRECISION_NO_GUARANTEE = -2.0;
	public static final double PRECISION_UNSPECIFIED = -1.0;

	private final T result;
	private final double precision;
	private final int numIters;
	private final double timeTaken;

	public ProcessingResult(T result, double timeTaken) {
		this(result, PRECISION_UNSPECIFIED, -1, timeTaken);
	}

	public ProcessingResult(T result, double precision, int numIters, double timeTaken)
	{
		assert timeTaken >= 0 && (precision >= 0 || precision == PRECISION_UNSPECIFIED || precision == PRECISION_NO_GUARANTEE);
		this.result = result;
		this.precision = precision;
		this.numIters = numIters;
		this.timeTaken = timeTaken;
	}

	public T getResult()
	{
		return result;
	}

	public double getPrecision()
	{
		return precision;
	}

	public int getNumIters()
	{
		return numIters;
	}

	public double getTimeTaken()
	{
		return timeTaken;
	}
}
