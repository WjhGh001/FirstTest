/**
 * @Date 2020-12-28
 * @author Junheng Wang
 */

public class ExaminationResults {
	private String xuehao;
	private String name;
	private String stuclass;
	private int grade;
	
	public ExaminationResults() {
		
	}
	public ExaminationResults(String xuehao, String name, String stuclass, int grade) {
		super();
		this.xuehao = xuehao;
		this.name = name;
		this.stuclass =stuclass;
		this.grade = grade;
	}
	
	public String getXuehao() {
		return xuehao;
	}
	public void setXuehao(String xuehao) {
		this.xuehao = xuehao;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	
	public String getStuclass() {
		return stuclass;
	}
	public void setStuclass(String stuclass) {
		this.stuclass = stuclass;
	}
	public int getGrade() {
		return grade;
	}
	public void setGrade(int grade) {
		this.grade = grade;
	}

	@Override
	public String toString() {
		return "WokerInfo [姓名 : " + name + ", 学号 : " + xuehao + ", 成绩 : " + grade + "]";
	}

}
