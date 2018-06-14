# Group-5
## **小組成員：**
* 物理19 鄭力珊<br />
* 物理20 黃紹倫<br />
* 物理20 李展誼<br />
* 物理20 葉百祥<br />
* 物理18 戴智恆<br />

## **簡介**
隨機震盪意指在特定非線性系統中增加雜訊，反而使訊雜比提升的特殊現象，目前觀察到的系統有：電路訊號系統、生物神經系統等。此次專題以神經為研究題材，探討之所以神經系統上可以產生隨機震盪的原因。

神經的運作模式為將外界的調幅訊號轉換為調頻訊號進行運算，而神經是一個極端非線性的系統且具有四個變數，十分難進行定量化分析，因此我們以目前生物領域常用來描述生物系統的參數：刺激大小與刺激差異。
我們可以很容易理解為什麼要研究刺激大小，但為什麼要看刺激差異呢?可以想像一隻青蛙正在水裡游泳，刺激大小意指他感覺得出熱不熱，刺激差異指鍋子裡的溫度變化的快不快，溫水煮青蛙之所以會把青蛙煮熟變是因為變是因為刺激差異太小所以他感覺不出來，也就是說刺激變化的快慢神經產生的訊號不一樣。

首先我們將x軸定為刺激強度，y軸定為刺激差異，畫出神經在不同條件底下產生的相圖。紅色部分標示出系統發生像變的位置，黃色箭頭為特定時間點雜訊將系統狀態從A點移動到B點：
![image](https://github.com/ShihPingLai/Group-5/blob/master/introduction.png)

一般認為隨機震盪發生在的相變附近(即可產生訊號與無法產生訊號的交界處)，又雜訊可同時趕變刺激大小與刺激差異，以往生物界認為神經之所以會在加入雜訊後產生訊號較佳的情況肇因於雜訊可以增加刺激強度，然而從圖中可發現同樣的雜訊造成的擾動反而對「刺激差異」此一因素關係較密切(兩者較垂直，因此只要增加一點點就可以產生相變)。

## **研究方向**
1.	SR
2.	改變雜訊的強度
3.	改變雜訊的頻率
4.	(I+noise – f) 的訊雜比
5.	Different type of noise
6.	(thermo phy)
7.	不同類型的訊號 sin wave/ square wave

## **Demo 方式**
1.  HH model 加雜訊
2.  音樂加雜訊<br />
  播放加雜訊後的音樂檔<br />
  Music資料夾
3.  Sine wave in rippletank simulation <br />
  以亮暗呈現類似水波槽俯視視角的動畫，分別展示沒有雜訊以及有雜訊的動畫圖形。<br />
  Ripple Tank Simulation資料夾

## 參考書籍：
*  ***Neuronal Dynamics***<br />
 ***From single neurons to networks and models of cognition***<br />
 *By Wulfram Gerstner, Werner M. Kistler, Richard Naud and Liam Paninski*<br />
  <http://neuronaldynamics.epfl.ch/>
  
*  ***Theoretical Neuroscience***<br />
 ***Computational and Mathematical Modeling of Neural Systems***<br />
 *By Laurence F. Abbott and Peter Dayan*<br />
  <https://mitpress.mit.edu/books/theoretical-neuroscience>
  
## **使用說明：**
### 水波槽模擬：<br />
需先將Ripple Tank Simulation資料夾中的檔案下載下來，分別有
1.  無雜訊
2.  0.1倍強度高斯白雜訊
3.  1倍強度高斯白雜訊
4.  2倍強度高斯白雜訊
5.  10倍強度高斯白雜訊<br />

(以 *Python 3* 執行檔案，即可看到模擬水波槽的動畫。)<br />

在大約加上2倍強度的高斯白雜訊時，水波槽的動畫會如同真實的水波一般，產生立體感。<br />
以應用來說，像是海面的影像，距離較遠的話僅能看出亮暗變化，若加上雜訊後使其產生隨機共振，可增強其對於海面的波浪辨識度。

### 雜訊函數: <br />
1.首先是直接以高斯函數產生亂數
2.另一種方式產生白雜訊，並將其家在方波輸入訊號

