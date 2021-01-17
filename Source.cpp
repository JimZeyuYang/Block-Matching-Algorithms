#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <chrono>

int addrY(int posnH, int posnW);
int addrU(int posnH, int posnW);
int addrV(int posnH, int posnW);
void fullSearch(int referenceH, int referenceW, int& resultH, int& resultW, int blockH,
	int blockW, unsigned char currentFrame[], unsigned char previousFrame[]);
void hexagonSearch(int referenceH, int referenceW, int& resultH, int& resultW, int blockH,
	int blockW, unsigned char currentFrame[], unsigned char previousFrame[]);
void diamondSearch(int referenceH, int referenceW, int& resultH, int& resultW, int blockH,
	int blockW, unsigned char currentFrame[], unsigned char previousFrame[]);
int SAD(int referenceH, int referenceW, int candidateH, int candidateW, int blockH, 
	int blockW, unsigned char currentFrame[], unsigned char previousFrame[]);
int SAD(int referenceH, int referenceW, int candidateH, int candidateW, int blockH,
	int blockW, int frameWidth, unsigned char currentFrame[], unsigned char previousFrame[]);
double calculatePSNR(unsigned char inputBuffer[], unsigned char outputBuffer[]);
int sobel(int referenceH, int referenceW, int blockH, int blockW, unsigned char currentFrame[]);
void AMMEA(int referenceH, int referenceW, int& resultH, int& resultW, int blockH, int blockW,
	unsigned char currentFrame[], unsigned char previousFrame[], unsigned char previousFrame4[],
	unsigned char previousFrame16[], unsigned char currentFrame4[], unsigned char currentFrame16[]);
void downSampler(unsigned char input[], unsigned char output[], int factor);
void fullSearch(int referenceH, int referenceW, int windowPosnH, int windowPosnW, int& resultH, int& resultW, int blockH, int blockW,
	int frameHeight, int frameWidth, int searchWindow, unsigned char currentFrame[], unsigned char previousFrame[]);
void L2FS(int referenceH, int referenceW, int resultH[3], int resultW[3], int blockH, int blockW, int frameHeight,
	int frameWidth, int searchWindow, unsigned char currentFrame[], unsigned char previousFrame[]);

int width = 832;
int height = 480;
int frames = 30;
int blockSize = 16;

int searchWindow = 64;

int frameSize = width * height * 1.5;

int main() {

	std::string algorithms[4] = {"Full Search", "Diamond Search", "Hexagon Search", "AMMEA"};

	std::ifstream file_org("BQMall_832x480_60.YUV", std::ios::binary);
	std::ofstream file_FS("Full_Search.YUV", std::ios::binary);
	std::ofstream file_DS("Diamond_Search.YUV", std::ios::binary);
	std::ofstream file_HS("Hexagon_Search.YUV", std::ios::binary);
	std::ofstream file_AMMEA("AMMEA.YUV", std::ios::binary);

	if (!file_org.is_open() || !file_DS.is_open() || !file_FS.is_open() || !file_AMMEA.is_open()) {
		std::cout << "The file couldn't be opened" << std::endl;
		exit(EXIT_FAILURE);
	}

	unsigned char* previousFrame = new unsigned char[frameSize];
	unsigned char* currentFrame = new unsigned char[frameSize];
	unsigned char* outputBuffer = new unsigned char[frameSize];

	unsigned char* previousFrame4 = new unsigned char[frameSize];
	unsigned char* previousFrame16 = new unsigned char[frameSize];
	unsigned char* currentFrame4 = new unsigned char[frameSize / 4 / 1.5];
	unsigned char* currentFrame16 = new unsigned char[frameSize / 16 / 1.5];

	double* PSNR = new double[frames];
	double* time = new double[frames];

	double avg_PSNR[4] = {};
	double tot_time[4] = {};
	double avg_time[4] = {};

	int matchH;
	int matchW;

	for (int A = 0; A < 4; A++) {
		std::cout << std::endl << "///// " << algorithms[A] << " /////" << std::endl << std::endl;

		for (double i = 0; i < frames; i++) {

			auto start = std::chrono::steady_clock::now();

			file_org.seekg(i * frameSize);
			file_org.read(reinterpret_cast<char*>(previousFrame), frameSize);
			file_org.seekg((i + 1) * frameSize);
			file_org.read(reinterpret_cast<char*>(currentFrame), frameSize);

			if (A == 3) {
				downSampler(previousFrame, previousFrame4, 2);
				downSampler(previousFrame, previousFrame16, 4);
				downSampler(currentFrame, currentFrame4, 2);
				downSampler(currentFrame, currentFrame16, 4);
			}


			for (int j = 0; j < height / blockSize; j++) {
				for (int k = 0; k < width / blockSize; k++) {
					int blockH = blockSize;
					int blockW = blockSize;
					if (j == height / blockSize) {
						blockH = height % blockSize;
					}
					if (k == width / blockSize) {
						blockW = width % blockSize;
					}
					
					if (A == 0) {
						fullSearch(j * blockSize, k * blockSize, matchH, matchW, blockH, blockW, currentFrame, previousFrame);
					}
					else if (A == 1) {
						diamondSearch(j * blockSize, k * blockSize, matchH, matchW, blockH, blockW, currentFrame, previousFrame);
					}
					else if (A == 2) {
						hexagonSearch(j * blockSize, k * blockSize, matchH, matchW, blockH, blockW, currentFrame, previousFrame);
					}
					else if (A == 3) {
						AMMEA(j * blockSize, k * blockSize, matchH, matchW, blockH, blockW,
							currentFrame, previousFrame, previousFrame4, previousFrame16, currentFrame4, currentFrame16);
					}
					

					for (int l = 0; l < blockH; l++) {
						for (int m = 0; m < blockW; m++) {
							outputBuffer[addrY(j * blockSize + l, k * blockSize + m)]
								= previousFrame[addrY(matchH + l, matchW + m)];

							outputBuffer[addrU(j * blockSize + l, k * blockSize + m)]
								= previousFrame[addrU(j * blockSize + l, k * blockSize + m)];

							outputBuffer[addrV(j * blockSize + l, k * blockSize + m)]
								= previousFrame[addrV(j * blockSize + l, k * blockSize + m)];
						}
					}
				}
			}

			PSNR[static_cast<int>(i)] = calculatePSNR(currentFrame, outputBuffer);

			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			time[static_cast<int>(i)] = elapsed_seconds.count();

			std::cout << "Frame " << i + 1 << "  |  PSNR: " << PSNR[static_cast<int>(i)];
			std::cout << "  |  Elapsed time: " << time[static_cast<int>(i)] << " sec" << std::endl;

			if (A == 0) {
				file_FS.seekp(i * frameSize);
				file_FS.write(reinterpret_cast<char*>(outputBuffer), frameSize);
			}
			else if (A == 1) {
				file_DS.seekp(i * frameSize);
				file_DS.write(reinterpret_cast<char*>(outputBuffer), frameSize);
			}
			else if (A == 2) {
				file_HS.seekp(i * frameSize);
				file_HS.write(reinterpret_cast<char*>(outputBuffer), frameSize);
			}
			else if (A == 3) {
				file_AMMEA.seekp(i * frameSize);
				file_AMMEA.write(reinterpret_cast<char*>(outputBuffer), frameSize);
			}

		}

		for (int i = 0; i < frames; i++) {
			avg_PSNR[A] += PSNR[i];
			tot_time[A] += time[i];
		}

		avg_PSNR[A] /= frames;
		avg_time[A] = tot_time[A] / frames;
	}

	std::cout << std::endl << "///////////" << std::endl;
	std::cout << "//Summary//" << std::endl;
	std::cout << "///////////" << std::endl << std::endl;

	for (int A = 0; A < 4; A++) {
		std::cout << algorithms[A] << std::endl;
		std::cout << "	Average PSNR: " << avg_PSNR[A] << std::endl;
		std::cout << "	Total elapsed time: " << tot_time[A] << " sec" << std::endl;
		std::cout << "	Average elapsed time per frame: " << avg_time[A] << " sec" << std::endl;
		std::cout << std::endl;
	}

	file_org.close();
	file_FS.close();
	file_DS.close();
	file_HS.close();
	file_AMMEA.close();

	delete[] outputBuffer;
	delete[] previousFrame;
	delete[] currentFrame;

	delete[] PSNR;
	delete[] time;
	
	return 0;
}

int addrY(int posnH, int posnW) {
	return posnH * width + posnW;
}

int addrU(int posnH, int posnW) {
	return height * width + posnH * width / 4 + posnW / 2;
}

int addrV(int posnH, int posnW) {
	return height * width * 1.25 + posnH * width / 4 + posnW / 2;
}

void fullSearch(int referenceH, int referenceW, int &resultH, int &resultW, int blockH, 
				int blockW, unsigned char currentFrame[], unsigned char previousFrame[]) {

	int windowH = referenceH - (searchWindow - blockH) / 2;
	int windowHH = referenceH + (searchWindow + blockH) / 2 - blockH;
	int windowW = referenceW - (searchWindow - blockW) / 2;
	int windowWW = referenceW + (searchWindow + blockW) / 2 - blockW;

	int minScore = 99999999;

	for (int i = windowH; i < windowHH; i++) {
		for (int j = windowW; j < windowWW; j++) {
			if (i >= 0 && i < height && j >= 0 && j < width) {
				int score = SAD(referenceH, referenceW, i, j, blockH, blockW, currentFrame, previousFrame);
				if (score < minScore) {
					minScore = score;
					resultH = i;
					resultW = j;
				}
			}
		}
	}
}

void hexagonSearch(int referenceH, int referenceW, int& resultH, int& resultW, int blockH,
	int blockW, unsigned char currentFrame[], unsigned char previousFrame[]) {

	int LHP[14] = { 0,0, 0,2, 2,1, 2,-1, 0,-2, -2,-1, -2,1 };
	int SHP[10] = { 0,0, 0,1, 1,0, 0,-1, -1,0 };

	int windowH = referenceH - (searchWindow - blockH) / 2;
	int windowHH = referenceH + (searchWindow + blockH) / 2 - blockH;
	int windowW = referenceW - (searchWindow - blockW) / 2;
	int windowWW = referenceW + (searchWindow + blockW) / 2 - blockW;

	bool center = false;
	int minScore;
	int centerH = referenceH;
	int centerW = referenceW;

	while (!center) {
		minScore = 99999999;
		for (int i = 0; i < 7; i++) {
			int testH = centerH + LHP[2 * i];
			int testW = centerW + LHP[2 * i + 1];
			if (testH >= 0 && testH <= height && testH >= windowH && testH <= windowHH
				&& testW >= 0 && testW <= width && testW >= windowW && testW <= windowWW) {
				int score = SAD(referenceH, referenceW, testH, testW, blockH, blockW, currentFrame, previousFrame);
				if (score < minScore) {
					minScore = score;
					resultH = testH;
					resultW = testW;
					if (i == 0) {
						center = true;
					}
					else {
						center = false;
					}
				}
			}
		}
		if (!center) {
			centerH = resultH;
			centerW = resultW;
		}
	}

	for (int i = 0; i < 5; i++) {
		int testH = referenceH + LHP[2 * i];
		int testW = referenceH + LHP[2 * i + 1];
		if (testH > 0 && testH < height - blockH && testH > windowH && testH < windowHH
			&& testW > 0 && testW < width - blockW && testW > windowW && testW < windowWW) {
			int score = SAD(referenceH, referenceW, testH, testW, blockH, blockW, currentFrame, previousFrame);
			if (score < minScore) {
				minScore = score;
				resultH = testH;
				resultW = testW;
			}
		}
	}
}

void diamondSearch(int referenceH, int referenceW, int& resultH, int& resultW, int blockH,
	int blockW, unsigned char currentFrame[], unsigned char previousFrame[]) {

	int LHP[18] = { 0,0, 0,2, 0,-2, 2,0, -2,0, 1,1, -1,1, 1,-1, -1,-1 };
	int SHP[10] = { 0,0, 0,1, 1,0, 0,-1, -1,0 };

	int windowH = referenceH - (searchWindow - blockH) / 2;
	int windowHH = referenceH + (searchWindow + blockH) / 2 - blockH;
	int windowW = referenceW - (searchWindow - blockW) / 2;
	int windowWW = referenceW + (searchWindow + blockW) / 2 - blockW;

	bool center = false;
	int minScore;
	int centerH = referenceH;
	int centerW = referenceW;

	while (!center) {
		minScore = 99999999;
		for (int i = 0; i < 9; i++) {
			int testH = centerH + LHP[2 * i];
			int testW = centerW + LHP[2 * i + 1];
			if (testH >= 0 && testH <= height && testH >= windowH && testH <= windowHH
				&& testW >= 0 && testW <= width && testW >= windowW && testW <= windowWW) {
				int score = SAD(referenceH, referenceW, testH, testW, blockH, blockW, currentFrame, previousFrame);
				if (score < minScore) {
					minScore = score;
					resultH = testH;
					resultW = testW;
					if (i == 0) {
						center = true;
					}
					else {
						center = false;
					}
				}
			}
		}
		if (!center) {
			centerH = resultH;
			centerW = resultW;
		}
	}

	for (int i = 0; i < 5; i++) {
		int testH = referenceH + LHP[2 * i];
		int testW = referenceH + LHP[2 * i + 1];
		if (testH > 0 && testH < height - blockH && testH > windowH && testH < windowHH
			&& testW > 0 && testW < width - blockW && testW > windowW && testW < windowWW) {
			int score = SAD(referenceH, referenceW, testH, testW, blockH, blockW, currentFrame, previousFrame);
			if (score < minScore) {
				minScore = score;
				resultH = testH;
				resultW = testW;
			}
		}
	}
}

int SAD(int referenceH, int referenceW, int candidateH, int candidateW, int blockH, 
				int blockW, unsigned char currentFrame[], unsigned char previousFrame[]) {
	int sum = 0;
	for (int i = 0; i < blockH; i++) {
		for (int j = 0; j < blockW; j++) {
			sum += abs(currentFrame[addrY(referenceH + i, referenceW + j)]
				- previousFrame[addrY(candidateH + i, candidateW + j)]);
		}
	}
	return sum;
}

int SAD(int referenceH, int referenceW, int candidateH, int candidateW, int blockH,
				int blockW, int frameWidth, unsigned char currentFrame[], unsigned char previousFrame[]) {
	int sum = 0;
	for (int i = 0; i < blockH; i++) {
		for (int j = 0; j < blockW; j++) {
			sum += abs(currentFrame[(referenceH + i) * frameWidth + referenceW + j]
				- previousFrame[(candidateH + i) * frameWidth + candidateW + j]);
		}
	}
	return sum;
}

double calculatePSNR(unsigned char inputBuffer[], unsigned char outputBuffer[]) {
	double MSE = 0;
	for (int i = 0; i <= height; i++) {
		for (int j = 0; j <= width; j++) {
			MSE += pow((inputBuffer[addrY(i, j)] - outputBuffer[addrY(i, j)]), 2);
		}
	}
	return 10 * log10(pow(155, 2) / (MSE / (static_cast<double>(width) * height)));
}

int sobel(int referenceH, int referenceW, int blockH, int blockW, unsigned char currentFrame[]) {
	int score = 0;
	int dx;
	int dy;
	for (int i = 1; i < blockH - 1; i++) {
		for (int j = 1; j < blockW - 1; j++) {
			dx = currentFrame[addrY(referenceH + i - 1, referenceW + j + 1)]
				+ 2 * currentFrame[addrY(referenceH + i, referenceW + j + 1)]
				+ currentFrame[addrY(referenceH + i + 1, referenceW + j + 1)]
				- currentFrame[addrY(referenceH + i - 1, referenceW + j - 1)]
				- 2 * currentFrame[addrY(referenceH + i, referenceW + j - 1)]
				- currentFrame[addrY(referenceH + i + 1, referenceW + j - 1)];

			dy = currentFrame[addrY(referenceH + i + 1, referenceW + j - 1)]
				+ 2 * currentFrame[addrY(referenceH + i + 1, referenceW + j)]
				+ currentFrame[addrY(referenceH + i + 1, referenceW + j + 1)]
				- currentFrame[addrY(referenceH + i - 1, referenceW + j - 1)]
				- 2 * currentFrame[addrY(referenceH + i - 1, referenceW + j)]
				- currentFrame[addrY(referenceH + i - 1, referenceW + j + 1)];

			score += abs(dx) + abs(dy);
		}
	}

	return score;
}

void AMMEA(int referenceH, int referenceW, int& resultH, int& resultW, int blockH, int blockW,
	unsigned char currentFrame[], unsigned char previousFrame[], unsigned char previousFrame4[],
	unsigned char previousFrame16[], unsigned char currentFrame4[], unsigned char currentFrame16[]) {

	int sad = SAD(referenceH, referenceW, referenceH, referenceW, blockH, blockW, width, currentFrame, previousFrame);

	if (sad < 500) {
		//LEVEL 0
		fullSearch(referenceH, referenceW, referenceH, referenceW, resultH, resultW, blockH, blockW, height, width, 24, currentFrame, previousFrame);
	}
	else {
		int homo = sobel(referenceH, referenceW, blockH, blockW, currentFrame);

		if (homo >= 20000) {
			//LEVEL 0
			fullSearch(referenceH, referenceW, referenceH, referenceW, resultH, resultW, blockH, blockW, height, width, 24, currentFrame, previousFrame);
		}
		else if (16000 <= homo && homo < 20000) {
			//LEVEL 1

			int H, W;
			fullSearch(referenceH / 2, referenceW / 2, referenceH / 2, referenceW / 2, H, W, blockH / 2, blockW / 2, height / 2, width / 2, 32, currentFrame4, previousFrame4);
			fullSearch(referenceH, referenceW, H * 2, W * 2, resultH, resultW, blockH, blockW, height, width, 16, currentFrame, previousFrame);

		}
		else if (homo < 16000) {
			//LEVEL 2

			int* H = new int[3];
			int* W = new int[3];
			int min = 99999999;
			int minH;
			int minW;
			L2FS(referenceH / 4, referenceW / 4, H, W, blockH / 4, blockW / 4, height / 4, width / 4, 30, currentFrame16, previousFrame16);

			for (int i = 0; i < 3; i++) {
				fullSearch(referenceH / 2, referenceW / 2, H[i] * 2, W[i] * 2, H[i], W[i], blockH / 2, blockW / 2, height / 2, width / 2, 16, currentFrame4, previousFrame4);
				int sad = SAD(referenceH / 2, referenceW / 2, H[i], W[i], blockH / 2, blockW / 2, width / 2, currentFrame4, previousFrame4);
				if (sad < min) {
					minH = H[i];
					minW = W[i];
				}
			}
			fullSearch(referenceH, referenceW, minH * 2, minW * 2, resultH, resultW, blockH, blockW, height, width, 8, currentFrame, previousFrame);
		}
	}
}

void downSampler(unsigned char input[], unsigned char output[], int factor) {
	for (int i = 0; i < height / factor; i++) {
		for (int j = 0; j < width / factor; j++) {
			output[i * width / factor + j] = input[i * width * factor + j * factor];
		}
	}

}

void fullSearch(int referenceH, int referenceW, int windowPosnH, int windowPosnW, int& resultH, int& resultW, int blockH, int blockW,
	int frameHeight, int frameWidth, int searchWindow, unsigned char currentFrame[], unsigned char previousFrame[]) {

	int windowH = windowPosnH - searchWindow / 2;
	int windowHH = windowPosnH + searchWindow / 2;
	int windowW = windowPosnW - searchWindow / 2;
	int windowWW = windowPosnW + searchWindow / 2;

	int minScore = 99999999;

	for (int i = windowH; i < windowHH; i++) {
		for (int j = windowW; j < windowWW; j++) {
			if (i >= 0 && i < frameHeight && j >= 0 && j < frameWidth) {
				int score = SAD(referenceH, referenceW, i, j, blockH, blockW, frameWidth, currentFrame, previousFrame);
				if (score < minScore) {
					minScore = score;
					resultH = i;
					resultW = j;
				}
			}
		}
	}
}

void L2FS(int referenceH, int referenceW, int resultH[3], int resultW[3], int blockH, int blockW, int frameHeight,
	int frameWidth, int searchWindow, unsigned char currentFrame[], unsigned char previousFrame[]) {

	int windowH = referenceH - (searchWindow - blockH) / 2;
	int windowHH = referenceH + (searchWindow + blockH) / 2 - blockH;
	int windowW = referenceW - (searchWindow - blockW) / 2;
	int windowWW = referenceW + (searchWindow + blockW) / 2 - blockW;

	int* minScore = new int[3];
	for (int i = 0; i < 3; i++) {
		minScore[i] = 99999999;
	}

	for (int i = windowH; i < windowHH; i++) {
		for (int j = windowW; j < windowWW; j++) {
			if (i >= 0 && i < frameHeight && j >= 0 && j < frameWidth) {
				int score = SAD(referenceH, referenceW, i, j, blockH, blockW, frameWidth, currentFrame, previousFrame);
				if (minScore[2] > score&& score >= minScore[1]) {
					minScore[2] = score;
					resultH[2] = i;
					resultW[2] = j;
				}
				else if (minScore[1] > score&& score >= minScore[0]) {
					minScore[2] = minScore[1];
					resultH[2] = resultH[1];
					resultW[2] = resultW[1];
					minScore[1] = score;
					resultH[1] = i;
					resultW[1] = j;
				}
				else if (minScore[0] > score) {
					minScore[2] = minScore[1];
					resultH[2] = resultH[1];
					resultW[2] = resultW[1];
					minScore[1] = minScore[0];
					resultH[1] = resultH[0];
					resultW[1] = resultW[0];
					minScore[0] = score;
					resultH[0] = i;
					resultW[0] = j;
				}
			}
		}
	}
}