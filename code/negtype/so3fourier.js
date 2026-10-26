// 길이 2 인 구조들은 회전 몫으로 SO(3) 의 점(상대 회전)이 되고, 둘 사이 D 는
// 두 상대 회전의 차이 각 θ 만의 함수 f(θ) = h·min(√6, 4 sin(θ/4)) 다 (h = 변 배율).
// 양쪽 불변 함수 f (f(e) = 0) 가 SO(3) 위에서 음의 유형 ⇔ 모든 l ≥ 1 에서
// 푸리에 계수 f̂_l = ∫ f χ_l dHaar ≤ 0. χ_l(θ) = sin((l+½)θ)/sin(θ/2), Haar 밀도 (1−cosθ)/π.
function coef(f, l, n = 200000) {
  let s = 0;
  for (let k = 0; k < n; k++) {
    const th = (k + 0.5) * Math.PI / n;
    const chi = Math.sin((l + 0.5) * th) / Math.sin(th / 2);
    s += f(th) * chi * (1 - Math.cos(th)) / Math.PI;
  }
  return s * Math.PI / n;
}
const fD = th => Math.min(Math.sqrt(6), 4 * Math.sin(th / 4));      // 실제 D (빈칸 상한 포함)
const fNoCap = th => 4 * Math.sin(th / 4);                            // 상한 없는 판 = 2·min(‖q−q′‖, ‖q+q′‖)
const fChord = th => Math.sqrt(4 * (1 - Math.cos(th)));            // 비교: 회전 행렬 현 거리 ‖R−R′‖_F (유클리드라 음의 유형)
for (const [name, f] of [['D (상한 포함)', fD], ['상한 없음', fNoCap], ['회전 행렬 현 거리', fChord]]) {
  const out = [];
  for (let l = 1; l <= 24; l++) out.push(coef(f, l));
  const pos = out.map((v, i) => [i + 1, v]).filter(([, v]) => v > 1e-9);
  console.log(name + ': 양수 계수 l =', pos.length ? pos.map(([l, v]) => l + '(' + v.toExponential(2) + ')').join(' ') : '없음');
  console.log('   f̂_1..8 =', out.slice(0, 8).map(v => v.toExponential(2)).join(' '));
}
