package simulator;

/**
 * Parameters for unbouded simulation methods.
 *
 * @author Przemysław Daca
 */
public class UnboundedSimulationParameters
{
	private UnboundedSimulationMethod unboudedSimulationMethod;
	private boolean simWhitebox;
	private double simFalseCnd;
	private double simMinProb;
	private double simRatio;
	private double terminationProb;
	private int simCheckBound;
	private int engine;

	public UnboundedSimulationParameters(UnboundedSimulationMethod usm, boolean simWhitebox, double simFalseCnd, double simMinProb, double simRatio,
			double terminationProb, int simCheckBound, int engine)
	{
		this.unboudedSimulationMethod = usm;
		this.simWhitebox = simWhitebox;
		this.simFalseCnd = simFalseCnd;
		this.simMinProb = simMinProb;
		this.simRatio = simRatio;
		this.terminationProb = terminationProb;
		this.simCheckBound = simCheckBound;
		this.engine = engine;
	}

	public UnboundedSimulationMethod getUnboudedSimulationMethod()
	{
		return unboudedSimulationMethod;
	}

	public boolean isSimWhitebox()
	{
		return simWhitebox;
	}

	public double getSimFalseCnd()
	{
		return simFalseCnd;
	}

	public double getSimMinProb()
	{
		return simMinProb;
	}

	public double getSimRatio()
	{
		return simRatio;
	}

	public double getTerminationProb()
	{
		return terminationProb;
	}

	public int getSimCheckBound()
	{
		return simCheckBound;
	}

	public int getEngine() {return engine; }

}
