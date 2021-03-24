package common;

import java.util.concurrent.TimeUnit;

public class Time
{
	private long start;

	public Time()
	{
		restart();
	}

	public void restart()
	{
		start = System.nanoTime();
	}

	public double elapsedSeconds(boolean restart) {
		if (!restart) {
			return elapsedSeconds();
		}
		double value = elapsedSeconds();
		restart();
		return value;
	}

	public double elapsedSeconds() {
		return (System.nanoTime() - start) / ((double) TimeUnit.SECONDS.toNanos(1));
	}
}
