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

package prism;

public class PrismDevNullLog extends PrismLog
{
	public PrismDevNullLog()
	{
	}

	@Override public boolean ready()
	{
		return true;
	}

	@Override public long getFilePointer()
	{
		return 0;
	}

	@Override public void flush()
	{
	}

	@Override public void close()
	{
	}

	// Basic print methods

	@Override public void print(boolean b)
	{
	}

	@Override public void print(char c)
	{
	}

	@Override public void print(double d)
	{
	}

	@Override public void print(float f)
	{
	}

	@Override public void print(int i)
	{
	}

	@Override public void print(long l)
	{
	}

	@Override public void print(Object obj)
	{
	}

	@Override public void print(String s)
	{
	}

	@Override public void println()
	{
	}
}

//------------------------------------------------------------------------------
