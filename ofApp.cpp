#include "ofApp.h"

void RungeKutta(MatrixXd &X, MatrixXd &dX, MatrixXd u, double dt, MatrixXd A, MatrixXd B) {
	MatrixXd k1 = A*X + B*u;
	MatrixXd k2 = A*(X + 0.5*k1*dt) + B*u;
	MatrixXd k3 = A*(X + 0.5*k2*dt) + B*u;
	MatrixXd k4 = A*(X + k3*dt) + B*u;
	dX = (k1 + 2.0*k2 + 2.0*k3 + k4)*dt / 6.0;
	X  += dX;
}

string serialize(const MatrixXd& M) {
	stringstream strm;

	for (int i = 0; i<M.rows(); i++) {
		for (int j = 0; j<M.cols(); j++) {
			strm << M(i, j);
			if (!(i == M.rows() - 1 && j == M.cols() - 1))
				strm << ",";
		}
		strm << "\n";
	}

	return strm.str();
}
double Uniform(void) {
	return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
}
double rand_normal(double mu, double sigma) {
	double z = sqrt(-2.0*log(Uniform())) * sin(2.0*M_PI*Uniform());
	return mu + sigma*z;
}


//--------------------------------------------------------------
void ofApp::setup(){
	ofSetVerticalSync(true); // sync with vertical refresh rate
	ofxGuiSetFont(ofToDataPath("ofxGraph/DIN Alternate Bold.ttf"), 10);

	graph.setup(100, 100, 800, 400);
	graph.setName("spring-mass-damper system wave");     // it automatically loads setting file, (sample.xml)
	graph.setDx(2.0); // which means delta of time
	graph.setColor(ofColor::white);  // ofColor(255,255,255)
	graph.setBufSize(1000);  
	graph.setLabel({ "y0","u" });

	frame = 0;

	gui.setup("parametor");
	gui.add(kk.setup("K", 1.0, 0.0, 5.0));
	gui.add(mm.setup("M", 0.1, 0.0, 1.0));
	gui.add(cc.setup("C", 0.1, 0.0, 1.0));

	AA.resize(2, 2);
	AA << 0, 1, -kk / mm, -cc / mm;

	BB.resize(2, 1);
	BB << 1, 1 / mm;

	CC.resize(1, 2);
	CC << 0, 1;

	DD.resize(1, 1);
	DD << 1;

	XX.resize(2, 1);
	XX(0, 0) = 0;
	XX(1, 0) = 0;
	dXX.resize(2, 1);
	dXX(0, 0) = 0;
	dXX(1, 0) = 0;

	uu.resize(1, 1);
	uu << 0;
	YY.resize(1, 1);
	YY(0, 0) = 0;


	dt = 0.01;
	tt = 0.0;



	Q_k.resize(2, 2);
	R_k.resize(1, 1);
	P_k.resize(2, 2);
	x_k_hat.resize(2, 1);
	x_k_hat_new.resize(2, 1);

	A_k = AA;
	C_k = CC;
	Q_k << .05, .05, .05, .05 ;
	R_k << 100.0;
	P_k << .0, 0, .0, 0;
	x_k_hat << 0, 0;
	x_k_hat_new << 0, 0;



}

//--------------------------------------------------------------
void ofApp::update() {
	frame++;
	tt += dt;

	/////////////////////
	AA << 0, 1, -kk / mm, -cc / mm;
	BB << 1, 1 / mm;

	double TT = 1000.0;
//	uu(0, 0) += sin(2.0*M_PI / TT*frame);
	//uu(0, 0) = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0)-0.5;


	RungeKutta(XX,dXX, uu, dt, AA, BB);
//	DD(0,0) = rand_normal(0, 1.0);
	YY = CC*XX+ DD;
	/////////////////////



//	kalman_filter
	x_k_hat_new = A_k * x_k_hat;
	P_k = A_k*P_k*A_k.transpose() + Q_k;
	K = P_k*C_k.transpose()*(C_k*P_k*C_k.transpose() + R_k).inverse();
	x_k_hat_new += K * (YY - C_k*x_k_hat_new);
	MatrixXd I = MatrixXd::Identity(A_k.rows(), A_k.rows());
	P_k = (I - K*C_k)*P_k;
	x_k_hat = x_k_hat_new;


	///////////////////////

	vector<float>value;
	value.push_back(YY(0, 0));
	//value.push_back(XX(1, 0));
	value.push_back(x_k_hat_new(1, 0));
	//value.push_back(XX(1, 0));
	//value.push_back(uu(0,0));

	graph.add(value);
	//graph2.add(-test);

}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackground(50, 50, 50);
	gui.draw();

	graph.draw();
	//graph2.draw();

	string info = "";
	info += "A mat = \n";
	info += serialize(AA);
	info += "B mat = \n";
	info += serialize(BB);
	info += "C mat = \n";
	info += serialize(CC);
	info += "D mat = \n";
	info += serialize(DD);

	string info2 = "";
	info2 += "X mat = \n";
	info2 += serialize(XX);
	info2 += "dX mat = \n";
	info2 += serialize(dXX);
	info2 += "u mat = \n";
	info2 += serialize(uu);
	info2 += "Y mat = \n";
	info2 += serialize(YY);


	string info3 = "";
	info3 += "x_hat mat = \n";
	info3 += serialize(x_k_hat);
	//info3 += "x_hat_new mat = \n";
	//info3 += serialize(x_k_hat_new);
	info3 += "P_k mat = \n";
	info3 += serialize(P_k);
	info3 += "K mat = \n";
	info3 += serialize(K);



	ofSetColor(240, 240, 240);
	ofFill();
	ofRect(20, ofGetHeight() * 2 /3, ofGetWidth()/2.0, ofGetHeight() / 3 -20);
	
	ofSetColor(240, 240, 0);
	if(uu(0, 0)==1)	ofSetColor(240, 0, 0);
	ofRect(ofGetWidth() / 2.0 + 60, ofGetHeight() * 2 / 3, ofGetWidth() / 3.0, ofGetHeight() / 3 - 20);

	ofSetColor(0);
	ofDrawBitmapString(info, 50, ofGetHeight() * 2 / 3 +40);
	ofDrawBitmapString(info2, ofGetWidth()/ 6.0, ofGetHeight() * 2 / 3 + 40);
	ofDrawBitmapString(info3, ofGetWidth()*2/ 6.0, ofGetHeight() * 2 / 3 + 40);

	string info_c="Control input";
	ofDrawBitmapString(info_c, ofGetWidth() / 2.0 + 60, ofGetHeight() * 2 / 3 + 40);


	ofSetColor(255);


}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == 's') {
		graph.saveSettings(); // save setting graph size and position
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
	uu(0, 0) = 1;
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
	uu(0, 0) = 0;

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
