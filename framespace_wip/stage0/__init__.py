"""0단계 코드 (연구 기획 E0 의 도구). 확정 패키지 framespace 위에 놓이고, E0 의 기준을 통과한 것만 패키지로 옮긴다.

- w2: 앙상블 W₂ (지상 비용 full / mean / band, 정확한 최적 수송), 바닥값, 편향 b/N, 부트스트랩, 삼각 부등식 점검, 방식 B 의 확률 대응
- decompose: M2′·M2″ 방식 A (맞춘 좌표의 평균·퍼짐 분해, 열과 변마다의 몫, 짝 없는 변의 정적 비용)
- band: 띠 제한 DP, 평균 점 출발
- sampling: 거리 자기상관과 간격 표본 (M1)
- esa: 탄성 모양 분석의 재구현 (비교 기준선)
- rulers: 자의 묶음 D_{w,ℓ} (위치 항) 와 고정 d₀ 의 TM 꼴 거리 D_TM, H-U2 의 요인 꼴, 삼각 부등식 위반률
"""
from .band import center_path, band_windows, dp_score_band, align_band, self_rotations, MeanStart
from .w2 import emd, w2_from_costs, scaled, pair_costs, ensemble_w2, soft_correspondence, replica_floors, half_split_floor, fit_bias, bootstrap_w2, triangle_violation
from .decompose import aligned_coordinates, decompose
from .sampling import paired_fixed_sq, distance_autocorr, spaced_indices, pooled_sample
from .esa import backbone_curve, resample, srvf, esa_distance, dp_reparam
from .rulers import edge_positions, d0_tm, start_poses, sub_fragment_poses, tm_align, tm_cost, wl_align, wl_cost, tm_forms, triangle_rate
